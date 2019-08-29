Name:  Casey Grove
Email: caseg1@umbc.edu
File:  seamcarving.cpp

Description: This program takes in an image, calculates the energy
	      (importance) of each pixel, and remove horizontal or
	      vertical rows of pixels based on the given output width
	      and height.

	     Problems with horizontal carving - slanted image is produced
	      when trying to shift pixels up. Doesn't find lowest energy
	      seam but rather starts with first lowest energy pixel and
	      creates seam from there.
	     
	     Program usage is:     
	      ./carve input_image output_width output_height
             
Help: http://www.cplusplus.com/reference/utility/pair/pair/
      https://stackoverflow.com/questions/7897050/adding-to-a-vector-of-pair
      https://stackoverflow.com/questions/9424173/find-the-smallest-amongst-3-numbers-in-c
      https://stackoverflow.com/questions/392932/how-do-i-use-the-conditional-operator
