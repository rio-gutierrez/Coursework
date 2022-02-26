#!/bin/bash

#this scripts takes the png files produced by gnuplot and turns them
# into an mp4 file

ffmpeg  -f image2  -i %07d.png -vcodec libx264 -crf 20  -pix_fmt yuv420p movie.mp4
