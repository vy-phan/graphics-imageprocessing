#include "convolution.h"

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <functional> // std::function
#include <vector>

using namespace graphics101;

namespace {
// Returns `value` clamped to lie in the range [min,max].
template< typename T >
inline T clamp( const T& value, const T& min, const T& max ) {
    return std::min( max, std::max( min, value ) );
}

// Your helper functions go here.

// NOTE: These are only suggestions!
// void convolve1D( int length, const ColorRGBA8* input, int input_stride, ColorRGBA8* output, int output_stride, const std::vector< int >& filter, bool normalize = true ) {}
// void resample1D( const ColorRGBA8* input, int input_length, int input_stride, ColorRGBA8* output, int output_length, int output_stride, const std::function<real(real)>& filter, real filter_radius, bool normalize = true ) {}

// To use this as a function taking only the x parameter, wrap it in a "lambda".
//     real radius = 2.0;
//     auto triangle_without_radius_parameter = [radius]( real x ) { return triangle( radius, x ); }
// Now we can call:
//     triangle_without_radius_parameter( x );
// or pass `triangle_without_radius_parameter` as an argument to a function taking a
// std::function<real(real)> parameter.

inline real triangle( real radius, real x ) {
    // Your code goes here.
    return std::max( 0., (1. - std::abs(x/radius)) );
}

}

namespace graphics101 {

// Converts each pixel to greyscale, saving the result into `output`.
void greyscale( const Image& input, Image& output ) {
    
    // Your code goes here.
    
    // Allocate an output image the same size as the input.
    output = Image( input.width(), input.height() );
    
    // Three ways to write it:
    // 1 Using input.pixel() and output.pixel()
    // 2 Using scanline() to iterate over rows.
    // 3 Using scanline() to iterate over columns.
    // I have written all three below.
    
    
    // 1 Using input.pixel()
    for( int j = 0; j < input.height(); ++j ) {
        for( int i = 0; i < input.width(); ++i ) {
            // The pixel() method gives direct access to the pixel at row i and column j.
            const ColorRGBA8 pix = input.pixel( i,j );
            
            // The grey value is just the average of red, green, and blue.
            // Perceptually, green is more important, but this is a fine approximation.
            const int grey = ( pix.r + pix.g + pix.b )/3;
            
            output.pixel( i,j ) = ColorRGBA8( grey, grey, grey );
        }
    }
    
    
    // 2 Using scanline() to iterate over each row.
    for( int j = 0; j < input.height(); ++j ) {
        
        // Get a pointer to the pixels in each row.
        // Pixels in a row are adjacent.
        // Pixels in a column are separated by a stride of input.width().
        const ColorRGBA8* row_in = input.scanline( j );
        const int stride_in = 1;
        
        ColorRGBA8* row_out = output.scanline( j );
        const int stride_out = 1;
        
        for( int i = 0; i < input.width(); ++i ) {
            // To access the i-th pixel, multiply i by the stride.
            const ColorRGBA8 pix = row_in[ i*stride_in ];
            
            // The grey value is just the average of red, green, and blue.
            // Perceptually, green is more important, but this is a fine approximation.
            const int grey = ( pix.r + pix.g + pix.b )/3;
            
            row_out[ i*stride_out ] = ColorRGBA8( grey, grey, grey );
        }
    }
    
    
    // 3 Using scanline() to iterate over each column.
    for( int i = 0; i < input.width(); ++i ) {
        
        // Get a pointer to the pixels in each column.
        // Pixels in a row are adjacent.
        // Pixels in a column are separated by a stride of input.width().
        const ColorRGBA8* col_in = input.scanline(0) + i;
        const int stride_in = input.width();
        
        ColorRGBA8* col_out = output.scanline(0) + i;
        const int stride_out = output.width();
        
        for( int j = 0; j < input.height(); ++j ) {
            // To access the j-th pixel, multiply j by the stride.
            const ColorRGBA8 pix = col_in[ j*stride_in ];
            
            // The grey value is just the average of red, green, and blue.
            // Perceptually, green is more important, but this is a fine approximation.
            const int grey = ( pix.r + pix.g + pix.b )/3;
            
            col_out[ j*stride_out ] = ColorRGBA8( grey, grey, grey );
        }
    }
}

// Subtracts `input1` from `input2`, saving the absolute value of the result into `output`.
// This function assumes that the dimensions of input1 and input2 match.
void difference( const Image& input1, const Image& input2, Image& output ) {
    assert( input1.width() == input2.width() );
    assert( input1.height() == input2.height() );
    
    // Your code goes here.
    
    // Allocate an output image the same size as the input.
    output = Image( input1.width(), input1.height() );
    
    // 1 Using input.pixel()
    for( int j = 0; j < input1.height(); ++j ) {
        for( int i = 0; i < input1.width(); ++i ) {
            const ColorRGBA8 pix1 = input1.pixel( i,j );
            const ColorRGBA8 pix2 = input2.pixel( i,j );
            
            const int rdiff = abs( int( pix1.r ) - int( pix2.r ) );
            const int gdiff = abs( int( pix1.g ) - int( pix2.g ) );
            const int bdiff = abs( int( pix1.b ) - int( pix2.b ) );
            
            output.pixel( i,j ) = ColorRGBA8( rdiff, gdiff, bdiff );
        }
    }
}

// Applies a box blur with `radius` to `input`, saving the result
// into `output`.
void blur_box( const Image& input, int radius, Image& output ) {
    assert( radius >= 0 );

    int w = input.width();
    int h = input.height();

    // Allocate an temp image the same size as the input.
    Image temp = Image( w, h );
    temp.fill( ColorRGBA8( 0, 0, 0 ) );

    // Allocate an output image the same size as the input.
    output = Image( w, h );

    // Vertically convolve to get an intermediate image
    for ( int x = 0; x < w; x++ ) {
        for( int y = 0; y < h; y++ ) {
            real count = 0.;
            real sumr = 0.;
            real sumg = 0.;
            real sumb = 0.;
            for ( int row = y - radius; row <= y + radius; row++ ) {
                if ( row < 0 || row >= h ) {
                    continue;
                }
                ColorRGBA8 pix = input.pixel( x, row );
                // red
                sumr = sumr + pix.r;
                // green
                sumg = sumg + pix.g;
                // blue
                sumb = sumb + pix.b;
                count++;
            }
            // write to intermediate image
            temp.pixel( x, y ) = ColorRGBA8( ( sumr / count ),
                                             ( sumg / count ),
                                             ( sumb / count ) );
        }
    }


    // Horizontally convolve with the intermediate image to get the output image
    for ( int x = 0; x < w; x++ ) {
        for ( int y = 0; y < h; y++ ) {
            real count = 0.;
            real sumr = 0.;
            real sumg = 0.;
            real sumb = 0.;
            for ( int column = x - radius; column <= x + radius; column++ ) {
                if ( column < 0 || column >= w ) {
                    continue;
                }
                ColorRGBA8 pix = temp.pixel( column, y );
                // red
                sumr = sumr + pix.r;
                // green
                sumg = sumg + pix.g;
                // blue
                sumb = sumb + pix.b;
                count++;
            }
            // write to output image
            output.pixel( x, y ) = ColorRGBA8( ( sumr/count ),
                                               ( sumg/count ),
                                               ( sumb/count ) );
        }
    }

}

// Scales the `input` image to the new dimensions, saving the result
// into `output`.
void scale( const Image& input, int new_width, int new_height, Image& output ) {
    assert( input.width() > 0 );
    assert( input.height() > 0 );
    assert( new_width > 0 );
    assert( new_height > 0 );

    int old_width = input.width();
    int old_height = input.height();

    // Allocate a temporary image that is scaled to the new width
    Image temp = Image( new_width, old_height );
    temp.fill( ColorRGBA8( 0, 0, 0 ) );

    // Allocate an output image with the new height and width
    output = Image( new_width, new_height );
    output.fill( ColorRGBA8( 0,0,0 ) );

    // Radius
    real xradius = real( old_width )  / real( new_width );
    real yradius = real( old_height ) / real( new_height );

    if ( new_width > old_width ) {
        xradius = 1.;
    }

    if ( new_height > old_height ) {
        yradius = 1.;
    }

    // Horizonally convolve to get an intermediate image
    for ( int y = 0; y < old_height; y++ ) {
        // each row of original image (old height)
        real xl = -0.5;
        real xh = old_width - 0.5;
        int n = new_width;
        real delta = ( xh - xl ) / real( n );
        real center = xl + ( delta / 2.0 );

        // Iterate through each pixel of the temp image with a new width
        for ( int i = 0; i <= n - 1; i++ ) {
            real sumr = 0.;
            real sumg = 0.;
            real sumb = 0.;
            real denominator = 0.;
            real x = center + i * delta;

            for ( int j = std::ceil( x - xradius ); j <= std::floor( x + xradius ); j++) {
                if ( j < 0 || j >= old_width ) {
                    continue;
                }

                ColorRGBA8 pix = input.pixel( j, y );
                real filter = triangle( xradius, x - j );
                // red
                sumr = sumr + pix.r * filter;
                // green
                sumg = sumg + pix.g * filter;
                // blue
                sumb = sumb + pix.b * filter;

                denominator += filter;
            }
            // Write to an intermediate image
            temp.pixel( i, y ) = ColorRGBA8( ( sumr / denominator ),
                                             ( sumg / denominator ),
                                             ( sumb / denominator ) );
        }
    }

    //Verically convolve from the temp image to get the output image
    for ( int x = 0; x < new_width; x++ ) {
        // each column of temporary image (new width)
        real yl = -0.5;
        real yh = old_height - 0.5;
        int n = new_height;
        real delta = ( yh - yl ) / real( n );
        real center = yl + ( delta / 2.0 );

        // Iterate through the pixels of the output image with a new height
        for ( int i = 0; i <= n - 1; i++ ) {
            real sumr = 0.;
            real sumg = 0.;
            real sumb = 0.;
            real denominator = 0.;
            real y = center + i * delta;

            for ( int j = std::ceil( y - yradius ); j <= std::floor( y + yradius ); j++ ) {
                if ( j < 0 || j >= old_height ) {
                    continue;
                }

                ColorRGBA8 pix = temp.pixel( x, j );
                real filter = triangle( yradius, y - j );
                // red
                sumr = sumr + pix.r * filter;
                // green
                sumg = sumg + pix.g * filter;
                // blue
                sumb = sumb + pix.b * filter;

                denominator += filter;
            }
            // Write to the output image
            output.pixel( x, i ) = ColorRGBA8( ( sumr / denominator ),
                                               ( sumg / denominator ),
                                               ( sumb / denominator ) );
        }
    }
}

// Convolves the `input` image with `filter`, saving the result into `output`.
// NOTE: This function assumes that `filter` is greyscale (has the same
//       values for red, green, and blue).
void convolve( const Image& input, const Image& filter, Image& output ) {
    // We assume that the filter is all gray.
    
    // Your code goes here.
    int rw = filter.width() / 2;
    int rh = filter.height() / 2;
    int w = input.width();
    int h = input.height();

    // Allocate an output image the same size as the input.
    output = Image( w, h );

    for ( int x = 0; x < w; x++ ) {
        for ( int y = 0; y < h; y++ ) {
            real sumr = 0.;
            real sumg = 0.;
            real sumb = 0.;
            real denominator = 0.;

            for ( int i = x - rw; i <= x + rw; i++ ) {
                for ( int j = y - rh; j <= y + rh; j++ ) {
                    if ( i < 0 || j < 0 || i >= w || j >= h ) {
                        continue;
                    }
                    const ColorRGBA8 pix = input.pixel(i, j);
                    const real value = filter.pixel( x - i + rw, y - j + rh ).r / 255.0;
                    // red
                    sumr = sumr + ( pix.r/255.0 ) * value;
                    // green
                    sumg = sumg + ( pix.g/255.0 ) * value;
                    // blue
                    sumb = sumb + ( pix.b/255.0 ) * value;

                    denominator += value;
                }
            }

            //normalize
            output.pixel( x, y ) = ColorRGBA8( ( sumr / denominator ) * 255.0,
                                               ( sumg / denominator ) * 255.0,
                                               ( sumb / denominator ) * 255.0 );
        }
    }
}

// Sharpens the `input` image by moving `amount` away from a blur with `radius`.
// Saves the result into `output`.
void sharpen( const Image& input, real amount, int radius, Image& output ) {
    assert( radius >= 0 );

    // Allocate an output image the same size as the input.
    output = Image( input.width(), input.height() );

    // Blur input image with box_blur
    Image blur_img;
    blur_box( input, radius, blur_img );

    // Your code goes here.
   for ( int x = 0; x < input.width(); x++ ) {
       for ( int y = 0; y < input.height(); y++ ) {
           ColorRGBA8 pix = input.pixel( x, y );
           ColorRGBA8 blr = blur_img.pixel( x, y );

           real r = ( 1. + amount ) * pix.r - amount * blr.r;
           real g = ( 1. + amount ) * pix.g - amount * blr.g;
           real b = ( 1. + amount ) * pix.b - amount * blr.b;
           // Write to output image
           output.pixel( x, y ) = ColorRGBA8( clamp( r, real(0.), real(255.0) ),
                                              clamp( g, real(0.), real(255.0) ),
                                              clamp( b, real(0.), real(255.0) ) );
       }
   }
}

// Performs edge detection on the `input` image. Stores the result into `output`.
void edge_detect( const Image& input, Image& output ) {

    int w = input.width();
    int h = input.height();

    // Allocate an output image the same size as the input.
    output = Image( w, h );

    // Edge Detection Filter:
    int filter[3] = { 1, 0, -1 };

    // Iterate through the image
    for ( int x = 0; x < w; x++ ) {
        for ( int y = 0; y < h; y++ ) {

            //Vertically 1D convolution
            int DyR = 0;
            int DyG = 0;
            int DyB = 0;
            for ( int j = y - 1; j <= y + 1; j++ ) {
                ColorRGBA8 pix;
                if ( j < 0 ) {
                    pix = input.pixel( x, j + 1 );
                }
                else if ( j >= h ) {
                    pix = input.pixel( x, j - 1 );
                }
                else {
                    pix = input.pixel( x, j );
                }

                DyR = DyR + pix.r * filter[ y - j + 1 ];
                DyG = DyG + pix.g * filter[ y - j + 1 ];
                DyB = DyB + pix.b * filter[ y - j + 1 ];
            }

            // Horizontally 1D convolution
            int DxR = 0;
            int DxG = 0;
            int DxB = 0;
            for ( int i = x - 1; i <= x + 1; i++ ) {
                ColorRGBA8 pix;
                if ( i < 0 ) {
                    pix = input.pixel( i + 1, y );
                }
                else if ( i >= w ) {
                    pix = input.pixel( i - 1, y );
                }
                else {
                    pix = input.pixel( i, y );
                }
                DxR = DxR + pix.r * filter[ x - i + 1 ];
                DxG = DxG + pix.g * filter[ x - i + 1 ];
                DxB = DxB + pix.b * filter[ x - i + 1 ];
            }

            // sqrt(Dx^2 + Dy^2)
            int r = int( std::round( std::sqrt( DxR * DxR + DyR * DyR ) ) );
            int g = int( std::round( std::sqrt( DxG * DxG + DyG * DyG ) ) );
            int b = int( std::round( std::sqrt( DxB * DxB + DyB * DyB ) ) );
            // Write to output image
            output.pixel( x, y ) = ColorRGBA8( clamp( r, 0, 255 ),
                                               clamp( g, 0, 255 ),
                                               clamp( b, 0, 255 ) );
        }
    }
}

}
