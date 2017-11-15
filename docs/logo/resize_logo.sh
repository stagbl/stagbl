#!/usr/bin/env bash
# Resize the logo and convert to PNG (for instance to use as a BitBucket avatar)
# requires ImageMagick
convert -resize 1400x1400 logo.pdf logo_resized.png
