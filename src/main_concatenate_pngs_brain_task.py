#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 18 17:54:44 2017

@author: ajoshi
"""

from PIL import Image


x_offset = 0
for ind in range(284):
    im1name = '/big_disk/ajoshi/coding_ground/brainsync/src/\
rest_vs_motor_after_rot_left_%d_m.png' % ind
    im2name = '/big_disk/ajoshi/coding_ground/brainsync/src/\
rest_vs_motor_after_rot_left_%d_d.png' % ind
    im3name = '/big_disk/ajoshi/coding_ground/brainsync/src/\
rest_vs_motor_after_rot_right_%d_m.png' % ind
    im4name = '/big_disk/ajoshi/coding_ground/brainsync/src/\
rest_vs_motor_after_rot_right_%d_d.png' % ind
    im5name = '/big_disk/ajoshi/coding_ground/brainsync/src/\
timing_%d.png' % ind

    im1 = Image.open(im1name)
    im2 = Image.open(im2name)
    im3 = Image.open(im3name)
    im4 = Image.open(im4name)
    im5 = Image.open(im5name)
    im5 = im5.resize((im1.size[0]+im2.size[0]+im3.size[0]+im4.size[0],
                      im5.size[1]), Image.ANTIALIAS)
    new_im = Image.new('RGB',
                       (im5.size[0], im1.size[1]+im5.size[1]))
    new_im.paste(im1, (0, 0))
    new_im.paste(im2, (im1.size[0], 0))
    new_im.paste(im3, (im1.size[0]+im2.size[0], 0))
    new_im.paste(im4, (im1.size[0]+im2.size[0]+im3.size[0], 0))
    new_im.paste(im5, (0, im1.size[1]))

    im1name = '/big_disk/ajoshi/coding_ground/brainsync/src/\
catimg_motor2rest_%d.png' % ind
    new_im.save(im1name)
