#include <random>
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

// Create a distribution
// uniform_real_distribution<double> distribution(0.0, 1.0);
// Use this to get a random value from our random engine e
// auto x = distribution(e);

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
    auto length = sqrt((x * x) + (y * y));
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

  for (unsigned int num_threads = 0; num_threads <= 6; ++num_threads) {
    auto total_threads = static_cast<unsigned int>(pow(2.0, num_threads));
    // Write number of threads
    cout << "Number of threads = " << total_threads << endl;
    // Write number of threads to the file
    data << "num_threads_" << total_threads;
    // Now execute 100 iterations
    for (unsigned int iters = 0; iters < 10; ++iters) {
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
      // We need to create total_threads threads
      vector<thread> threads;
      for (unsigned int n = 0; n < total_threads; ++n)
        // Working in base 2 to make things a bit easier
        threads.push_back(
            thread(monte_carlo_pi,
                   static_cast<unsigned int>(pow(2.0, 24.0 - num_threads))));
      // Join the threads (wait for them to finish)
      for (auto &t : threads)
        t.join();
      // Get the end time
      auto end = system_clock::now();
      // Get the total time
      auto total = end - start;
      // Convert to milliseconds and output to file
      double pi = (4.0 * total_circle) / static_cast<double>(total_points);
      auto ms = duration_cast<milliseconds>(total).count();
      double actualPI = 3.1415926535897931;
      auto error = abs((pi - actualPI) / actualPI) * 100.0;
      cout << total_threads << " threads computed Pi: " << pi << " In: " << ms
           << " milliseconds, Error: " << error << "%" << endl;
      data << ", " << ms;
    }
    data << endl;
  }
  // Close the file
  data.close();
  return 0;
}