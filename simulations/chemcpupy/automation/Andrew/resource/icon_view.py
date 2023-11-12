#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: chrisarcadia & jacobrosenstein
@created: 2018/03/22
"""

import io
import matplotlib
import base64

filename = "./resource/iconBase64Image.txt";

with open(filename, "rb") as file:
    imageBase64 = file.read();


def viewIcon(imageBase64):
    image = base64.decodebytes(imageBase64)
    imageIO = io.BytesIO(image)
    imageRead = matplotlib.image.imread(imageIO, format='JPG')  
    matplotlib.pyplot.imshow(imageRead, interpolation='nearest')
    matplotlib.pyplot.axis('off')
    matplotlib.pyplot.show()
    
viewIcon(imageBase64);
    