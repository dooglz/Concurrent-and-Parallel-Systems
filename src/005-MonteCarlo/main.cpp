#include <random>
#include <String>
#include <thread>
#include <chrono>
#include <iostream>
#include <fstream>
#include <mutex>

using namespace std;
using namespace std::chrono;

unsigned int total_points;
unsigned int total_circle;
mutex total_points_mutex;
mutex total_circle_mutex;

#define MIN_THREADS 8
#define MAX_THREADS 8
#define LOOPS 10
#define MIN_POWER 8
#define MAX_POWER 32
// Create a distribution
// uniform_real_distribution<double> distribution(0.0, 1.0);
// Use this to get a random value from our random engine e
// auto x = distribution(e);

template <typename T> T average(T t[], int n) {
  T s = t[n - 1];
  for (int i = 0; i < (n - 1); i++)
    s += t[i];
  return s / n;
}
unsigned int correct_pi_digits(double p) {
  double actualPI = 3.1415926535897931;
  unsigned int i;
  for (i = 1; i < 16; i++) {
    double multiple = pow(10.0, i);
    double a = trunc(multiple * p);
    double b = trunc(multiple * actualPI);
    if (a != b)
      break;
  }
  return i;
}

void monte_carlo_pi(unsigned int iterations) {
  // Create a random engine
  auto millis =
      duration_cast<milliseconds>(system_clock::now().time_since_epoch());
  default_random_engine e(millis.count());
  // Create a distribution - we want doubles between 0.0 and 1.0
  uniform_real_distribution<double> distribution(0.0, 1.0);

  // Keep track of number of points in circle
  unsigned int in_circle = 0;
  // Iterate
  for (unsigned int i = 0; i < iterations; ++i) {
    // Generate random point
    auto x = distribution(e);
    auto y = distribution(e);
    // Get length of vector defined - use Pythagarous
    auto length = abs(x * x) + abs(y * y);
    // Check if in circle
    if (length <= 1.0)
      ++in_circle;
  }
  {
    std::lock_guard<std::mutex> lock(total_points_mutex);
    total_points += iterations;
  }
  {
    std::lock_guard<std::mutex> lock(total_circle_mutex);
    total_circle += in_circle;
  }
  // Calculate pi
  // auto pi = (4.0 * in_circle) / static_cast<double>(iterations);
  // cout << pi << endl;
}

int main() {
  // Create data file
  ofstream data("montecarlo.csv", ofstream::out);
  data << "numbers per thread, threads, Avg Time(ms), Avg PI, Avg Error"
       << endl;
  for (unsigned int total_numbers = 1024; true;
       total_numbers+= 1024) {
    cout << "\n total_numbers = " << total_numbers << " - "
         << total_numbers << endl;
    for (unsigned int num_threads = MIN_THREADS; num_threads <= MAX_THREADS;
         ++num_threads) {
      unsigned int thread_iterations =
          static_cast<unsigned int>(total_numbers / num_threads);
      if (thread_iterations < 3) {
        continue;
      }
      // Write number of threads
      // cout << "Number of threads = " << num_threads << ", Numbers per thread
      // = " << thread_iterations << endl;
      // Write number of threads to the file
      unsigned int pi_digits[LOOPS];
      long long times[LOOPS];
      double errors[LOOPS];
      for (unsigned int iters = 0; iters < LOOPS; ++iters) {
        {
          std::lock_guard<std::mutex> lock(total_points_mutex);
          total_points = 0;
        }
        {
          std::lock_guard<std::mutex> lock(total_circle_mutex);
          total_circle = 0;
        }

        // Get the start time
        auto start = system_clock::now();

        //----------------------- Do work -----------------
        // We need to create total_threads threads
        vector<thread> threads;
        for (unsigned int n = 0; n < num_threads; ++n)
          // Working in base 2 to make things a bit easier
          threads.push_back(thread(monte_carlo_pi, thread_iterations));
        // Join the threads (wait for them to finish)
        for (auto &t : threads)
          t.join();
        //--------------------------------------------------

        // Get the end time
        auto end = system_clock::now();
        // Get the total time
        auto total = end - start;
        // Convert to milliseconds and output to file
        double pi = (4.0 * total_circle) / static_cast<double>(total_points);
        long long ms = duration_cast<milliseconds>(total).count();
        double actualPI = 3.1415926535897931;
        double error = abs((pi - actualPI) / actualPI) * 100.0;
        pi_digits[iters] = correct_pi_digits(pi);
        times[iters] = ms;
        errors[iters] = error;

        /*
        cout << num_threads << " threads computed Pi: " << pi << " In: " << ms
             << " milliseconds, Error: " << error << "%" << endl;
        data << ", " << ms << "ms";
        */
      }
      auto avg_digits = floor(average(pi_digits, LOOPS));
      auto avg_time = average(times, LOOPS);
      auto avg_err = average(errors, LOOPS);
      cout << num_threads << " threads\tnumbers " << thread_iterations << "\t"
           << avg_time << "ms\tscore " << avg_digits << " "
           << std::string(avg_digits, '#') << endl;
      data << total_numbers << "," << num_threads << "," << avg_time << ","
           << avg_digits << "," << avg_err << endl;
      data << endl;
    }
  }
  // Close the file
  data.close();
  return 0;
}