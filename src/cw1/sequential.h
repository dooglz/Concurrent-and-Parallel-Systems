extern const unsigned int SIZE;
extern const double CHECK_VALUE;

namespace seq {
  int start();
}
#include <fstream>
#include <iostream>
#include <string>

namespace intarray2bmp
{

  //-------------------------------------------------------------------------- 
  // This little helper is to write little-endian values to file.
  //
  struct lwrite
  {
    unsigned long value;
    unsigned      size;
    lwrite(unsigned long value, unsigned size) :
      value(value), size(size)
    { }
  };

  //--------------------------------------------------------------------------
  inline std::ostream& operator << (std::ostream& outs, const lwrite& v)
  {
    unsigned long value = v.value;
    for (unsigned cntr = 0; cntr < v.size; cntr++, value >>= 8)
      outs.put(static_cast <char> (value & 0xFF));
    return outs;
  }

  //--------------------------------------------------------------------------
  // Take an integer array and convert it into a color image.
  //
  // This first version takes an array of array style of array:
  //   int* a[ 10 ]
  //
  // The second, overloaded version takes a flat C-style array:
  //   int a[ 10 ][ 10 ]
  //
  template <typename doubleType>
  bool intarray2bmp(
    const std::string& filename,
    doubleType**          intarray,
    unsigned           rows,
    unsigned           columns,
    doubleType            min_value,
    doubleType            max_value
    ) {
    // This is the difference between each color based upon
    // the number of distinct values in the input array.
    double granularity = 360.0 / ((double)(max_value - min_value) + 1);

    // Open the output BMP file
    std::ofstream f(filename.c_str(),
      std::ios::out | std::ios::trunc | std::ios::binary);
    if (!f) return false;

    // Some basic
    unsigned long headers_size = 14  // sizeof( BITMAPFILEHEADER )
      + 40; // sizeof( BITMAPINFOHEADER )
    unsigned long padding_size = (4 - ((columns * 3) % 4)) % 4;
    unsigned long pixel_data_size = rows * ((columns * 3) + padding_size);

    // Write the BITMAPFILEHEADER
    f.put('B').put('M');                           // bfType
    f << lwrite(headers_size + pixel_data_size, 4);  // bfSize
    f << lwrite(0, 2);  // bfReserved1
    f << lwrite(0, 2);  // bfReserved2
    f << lwrite(headers_size, 4);  // bfOffBits

    // Write the BITMAPINFOHEADER
    f << lwrite(40, 4);  // biSize
    f << lwrite(columns, 4);  // biWidth
    f << lwrite(rows, 4);  // biHeight
    f << lwrite(1, 2);  // biPlanes
    f << lwrite(24, 2);  // biBitCount
    f << lwrite(0, 4);  // biCompression=BI_RGB
    f << lwrite(pixel_data_size, 4);  // biSizeImage
    f << lwrite(0, 4);  // biXPelsPerMeter
    f << lwrite(0, 4);  // biYPelsPerMeter
    f << lwrite(0, 4);  // biClrUsed
    f << lwrite(0, 4);  // biClrImportant

    // Write the pixel data
    for (unsigned row = rows; row; row--)           // bottom-to-top
    {
      for (unsigned col = 0; col < columns; col++)  // left-to-right
      {
      unsigned char value;
      double        d = 0.0;
      //d =  d - min_value * granularity;
    
      // If we haven't overrun the end of our input, convert it to a grayscale value.
      // Input is clamped to the range [0, 1], where 0 --> black and 1 --> white.
    //  cout << (row - 1) << " - " << col << endl;
      d = intarray[row - 1][col];
      d = (d + (max_value - (min_value + max_value))) / (max_value - min_value);
        if (d < 0.0) d = 0.0;
        else if (d > 1.0) d = 1.0;

        value = 255 * d;

      f.put(static_cast <char> (value))
       .put(static_cast <char> (value))
        .put(static_cast <char> (value));


    }
     // cout << endl;
    if (padding_size) f << lwrite(0, padding_size);
  }


    // All done!
    return f.good();
  }

  //--------------------------------------------------------------------------
 /* template <typename doubleType>
  bool intarray2bmp(
    const std::string& filename,
    doubleType*           intarray,
    unsigned           rows,
    unsigned           columns,
    doubleType            min_value,
    doubleType            max_value
    ) {
    doubleType** ia = new(std::nothrow) doubleType*[rows];
    for (unsigned row = 0; row < rows; row++)
    {
      ia[row] = intarray + (row * columns);
    }
    bool result = intarray2bmp(
      filename, ia, rows, columns, min_value, max_value
      );
    delete[] ia;
    return result;
  }
*/
} // namespace intarray2bmp