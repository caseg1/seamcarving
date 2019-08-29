#define cimg_display 0
#include "CImg.h"
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <utility>
using namespace cimg_library;


int main(int argc, char *argv[]) {
  CImg<double> input(argv[1]);
  CImg<double> lab = input.RGBtoLab();
  Eigen::Vector3d *image = new Eigen::Vector3d[input.width()*input.height()];
  for (unsigned int i=0; i<input.width(); i++) {
    for (unsigned int j=0; j<input.height(); j++) {
      image[i*input.height()+j][0] = lab(i, j, 0);
      image[i*input.height()+j][1] = lab(i, j, 1);
      image[i*input.height()+j][2] = lab(i, j, 2);
    }
  }
  
  float energy[input.width()*input.height()];  //array to hold pixel energy 
  float dx,dy;
  
  for (unsigned int i=0; i<input.width(); i++) {
    for (unsigned int j=0; j<input.height(); j++) {
      
      //corners
      if (i == 0 && j == 0) {//top left corner
	dx = (image[(1+i)*input.height()+j] - image[i*input.height()+j]).norm() / 1.0;
	dy = (image[i*input.height()+(j+1)] - image[i*input.height()+j]).norm() / 1.0;
      }
      else if (i == 0 && j == input.height()-1) {//bottom left corner
	dx = (image[(i+1)*input.height()+(j+input.height()-1)] - image[i*input.height()+j]).norm() / 1.0;
	dy = (image[i*input.height()+(j+input.height()-2)] - image[i*input.height()+j]).norm() / 1.0;
      }
      else if (i == input.width()-1 && j == 0) {//top right corner
	dx = (image[(input.width()-2+i)*input.height()+(j+input.height())] - image[i*input.height()+j]).norm() / 1.0;
	dy = (image[(input.width()-1+i)*input.height()+(j+input.height()+1)] - image[i*input.height()+j]).norm() / 1.0;
      }
      else if (i == input.width()-1 && j == input.height()-1) {//bottom right corner
	dx = (image[(input.width()-2+i)*input.height()+(j+input.height()-1)] - image[i*input.height()+j]).norm() / 1.0;
	dy = (image[(input.width()-1+i)*input.height()+(j+input.height()-2)] - image[i*input.height()+j]).norm() / 1.0;
      }
      
      //  edges
      else if ( i == 0 && (j > 0 && j < input.height()-1) ) {//left
	dx = (image[(i+1)*input.height()+j] - image[i*input.height()+j]).norm() / 1.0;
	dy = (image[i*input.height()+(j-1)] - image[i*input.height()+(j+1)]).norm() / 2.0;
      }
      else if ( i == input.width()-1 && (j > 0 && j < input.height()-1) ) {//right
	dx = (image[(i-1)*input.height()+j] - image[i*input.height()+j]).norm() / 1.0;
	dy = (image[i*input.height()+(j-1)] - image[i*input.height()+(j+1)]).norm() / 2.0;
      }
      else if ( j == 0 && (i > 0 && i < input.width()-1) ) {//top
	dx = (image[(i-1)*input.height()+j] - image[(i+1)*input.height()+j]).norm() / 2.0;
	dy = (image[i*input.height()+(j+1)] - image[i*input.height()+j]).norm() / 1.0;
      }
      else if ( j == 0 && (i > 0 && i < input.width()-1) ) {//bottom
	dx = (image[(i-1)*input.height()+j] - image[(i+1)*input.height()+j]).norm() / 2.0;
	dy = (image[i*input.height()+(j+1)] - image[i*input.height()+j]).norm() / 1.0;
      }
      else {//insides
	dx = (image[(i-1)*input.height()+j] - image[(i+1)*input.height()+j]).norm() / 2.0;
	dy = (image[i*input.height()+(j-1)] - image[i*input.height()+(j+1)]).norm() / 2.0;;
      }
      energy[i*input.height()+j] = (dx + dy);
    }
  }
  
  //Carve
  CImg<double> output(atoi(argv[3]), atoi(argv[4]), input.depth(), input.spectrum(), 0);

  //find number of vseams/hseams to remove
  int vseams = input.width() - output.width();
  int hseams = input.height() - output.height();

  while (vseams) {

    std::vector<std::pair<int,int>> seam;
    int x;
    float min_e = 10000.0;

    for (unsigned int i=0; i<input.width(); i++) {
      if (energy[i*input.height()] < min_e) {
	min_e = energy[i*input.height()+0];
	x = i;
      }
    }
    seam.push_back(std::make_pair(x,0));

    for (unsigned int j=1; j<input.height(); j++) {
      if (x == 0) {//left edge	
	energy[x*input.height()+j] < energy[(x+1)*input.height()+j] ? x : x++;
      }
      else if (x==input.width()-1) {//right edge
	energy[x*input.height()+j] < energy[(x-1)*input.height()+j] ? x : x--;
      }
      else {//non-edge pixel
	energy[x*input.height()+j] < energy[(x-1)*input.height()+j] ? x : x--;
	energy[x*input.height()+j] < energy[(x+1)*input.height()+j] ? x : x++;
      }
      seam.push_back(std::make_pair(x,j));

      //red line of seam for debugging
      image[x*input.height()+j][0] = 0.0;
      image[x*input.height()+j][1] = 255.0;
      image[x*input.height()+j][2] = 255.0;
    }   

    //loop through seam[(x,y)] x=?, y=0,1,2,3,4...
    for (unsigned int j=0; j<input.height(); j++) {
      for (unsigned int i=seam[j].first; i<input.width(); i++) {
	image[i*input.height()+j] = image[(i+1)*input.height()+j];
      }
    }    
    
    //redo energy function for nearby pixels
    for (unsigned int j=0; j<input.width(); j++) {
      for (unsigned int i=seam[j].first-1; i<seam[j].first+2; i++) {
	
	//corners
	if (i == 0 && j == 0) {//top left corner
	  dx = (image[(1+i)*input.height()+j] - image[i*input.height()+j]).norm() / 1.0;
	  dy = (image[i*input.height()+(j+1)] - image[i*input.height()+j]).norm() / 1.0;
	}
	
	else if (i == 0 && j == input.height()-1) {//bottom left corner
	  dx = (image[(i+1)*input.height()+(j+input.height()-1)] - image[i*input.height()+j]).norm() / 1.0;
	  dy = (image[i*input.height()+(j+input.height()-2)] - image[i*input.height()+j]).norm() / 1.0;
	}
	else if (i == input.width()-1 && j == 0) {//top right corner
	  dx = (image[(input.width()-2+i)*input.height()+(j+input.height())] - image[i*input.height()+j]).norm() / 1.0;
	  dy = (image[(input.width()-1+i)*input.height()+(j+input.height()+1)] - image[i*input.height()+j]).norm() / 1.0;
	}
	else if (i == input.width()-1 && j == input.height()-1) {//bottom right corner
	  dx = (image[(input.width()-2+i)*input.height()+(j+input.height()-1)] - image[i*input.height()+j]).norm() / 1.0;
	  dy = (image[(input.width()-1+i)*input.height()+(j+input.height()-2)] - image[i*input.height()+j]).norm() / 1.0;
	}
	
	//  edges
	else if ( i == 0 && (j > 0 && j < input.height()-1) ) {//left
	  dx = (image[(i+1)*input.height()+j] - image[i*input.height()+j]).norm() / 1.0;
	  dy = (image[i*input.height()+(j-1)] - image[i*input.height()+(j+1)]).norm() / 2.0;
	}
	else if ( i == input.width()-1 && (j > 0 && j < input.height()-1) ) {//right
	  dx = (image[(i-1)*input.height()+j] - image[i*input.height()+j]).norm() / 1.0;
	  dy = (image[i*input.height()+(j-1)] - image[i*input.height()+(j+1)]).norm() / 2.0;
	}
	else if ( j == 0 && (i > 0 && i < input.width()-1) ) {//top
	  dx = (image[(i-1)*input.height()+j] - image[(i+1)*input.height()+j]).norm() / 2.0;
	  dy = (image[i*input.height()+(j+1)] - image[i*input.height()+j]).norm() / 1.0;
	}
	else if ( j == 0 && (i > 0 && i < input.width()-1) ) {//bottom
	  dx = (image[(i-1)*input.height()+j] - image[(i+1)*input.height()+j]).norm() / 2.0;
	  dy = (image[i*input.height()+(j+1)] - image[i*input.height()+j]).norm() / 1.0;
	}
	else {//insides
	  dx = (image[(i-1)*input.height()+j] - image[(i+1)*input.height()+j]).norm() / 2.0;
	  dy = (image[i*input.height()+(j-1)] - image[i*input.height()+(j+1)]).norm() / 2.0;;
	}
	energy[i*input.height()+j] = (dx + dy);
      }
    }
    seam.clear();
    vseams--;
  }
  
  while (hseams) {

    std::vector<std::pair<int,int>> seam;
    int y;
    float min_e = 10000.0;
    
    //find first pixel
    for (unsigned int j=0; j<input.height(); j++) {
      if (energy[input.height()+j] < min_e) {
	min_e = energy[j*input.height()];
	y = j;
      }
    }
    seam.push_back(std::make_pair(0,y));
    
    for (unsigned int i=1; i<input.width(); i++) {
      if (y == 0) {//top edge
        energy[i*input.height()+y] < energy[(i+1)*input.height()+y] ? y : y++;
      }
      else if (y==input.width()-1) {//bottom edge
        energy[i*input.height()+y] < energy[(i-1)*input.height()+y] ? y : y--;
      }
      else {//non-edge pixel
        energy[i*input.height()+y] < energy[(i-1)*input.height()+y] ? y : y--;
        energy[i*input.height()+y] < energy[(i+1)*input.height()+y] ? y : y++;
      }
      seam.push_back(std::make_pair(i,y));

      //red line of seam for debugging
      image[i*input.height()+y][0] = 0.0;
      image[i*input.height()+y][1] = 255.0;
      image[i*input.height()+y][2] = 255.0;
    }

    //!! BUGGY !! produces slanted picture
    //loop through seam[(x,y)] x=0,1,2,3..., y=?
    for (unsigned int i=0; i<input.width(); i++) {
      for (unsigned int j=seam[i].second; j<input.height(); j++) {
	image[i*input.height()+j] = image[i*input.height()+(j+1)];
      }
    }
    
    //redo energy function for nearby pixels
    for (unsigned int i=0; i<input.height(); i++) {
      for (unsigned int j=seam[i].second-1; j<seam[i].second+2; j++) {
	
	//corners
	if (i == 0 && j == 0) {//top left corner
	  dx = (image[(1+i)*input.height()+j] - image[i*input.height()+j]).norm() / 1.0;
	  dy = (image[i*input.height()+(j+1)] - image[i*input.height()+j]).norm() / 1.0;
	}
	
	else if (i == 0 && j == input.height()-1) {//bottom left corner
	  dx = (image[(i+1)*input.height()+(j+input.height()-1)] - image[i*input.height()+j]).norm() / 1.0;
	  dy = (image[i*input.height()+(j+input.height()-2)] - image[i*input.height()+j]).norm() / 1.0;
	}
	else if (i == input.width()-1 && j == 0) {//top right corner
	  dx = (image[(input.width()-2+i)*input.height()+(j+input.height())] - image[i*input.height()+j]).norm() / 1.0;
	  dy = (image[(input.width()-1+i)*input.height()+(j+input.height()+1)] - image[i*input.height()+j]).norm() / 1.0;
	}
	else if (i == input.width()-1 && j == input.height()-1) {//bottom right corner
	  dx = (image[(input.width()-2+i)*input.height()+(j+input.height()-1)] - image[i*input.height()+j]).norm() / 1.0;
	  dy = (image[(input.width()-1+i)*input.height()+(j+input.height()-2)] - image[i*input.height()+j]).norm() / 1.0;
	}
	
	//  edges
	else if ( i == 0 && (j > 0 && j < input.height()-1) ) {//left
	  dx = (image[(i+1)*input.height()+j] - image[i*input.height()+j]).norm() / 1.0;
	  dy = (image[i*input.height()+(j-1)] - image[i*input.height()+(j+1)]).norm() / 2.0;
	}
	else if ( i == input.width()-1 && (j > 0 && j < input.height()-1) ) {//right
	  dx = (image[(i-1)*input.height()+j] - image[i*input.height()+j]).norm() / 1.0;
	  dy = (image[i*input.height()+(j-1)] - image[i*input.height()+(j+1)]).norm() / 2.0;
	}
	else if ( j == 0 && (i > 0 && i < input.width()-1) ) {//top
	  dx = (image[(i-1)*input.height()+j] - image[(i+1)*input.height()+j]).norm() / 2.0;
	  dy = (image[i*input.height()+(j+1)] - image[i*input.height()+j]).norm() / 1.0;
	}
	else if ( j == 0 && (i > 0 && i < input.width()-1) ) {//bottom
	  dx = (image[(i-1)*input.height()+j] - image[(i+1)*input.height()+j]).norm() / 2.0;
	  dy = (image[i*input.height()+(j+1)] - image[i*input.height()+j]).norm() / 1.0;
	}
	else {//insides
	  dx = (image[(i-1)*input.height()+j] - image[(i+1)*input.height()+j]).norm() / 2.0;
	  dy = (image[i*input.height()+(j-1)] - image[i*input.height()+(j+1)]).norm() / 2.0;;
	}
	energy[i*input.height()+j] = (dx + dy);
      }
    }
    seam.clear();
    hseams--;
  }
  
  float m = 0.0;
  for (unsigned int i=0; i<output.width(); i++) {
    for (unsigned int j=0; j<output.height(); j++) {
      m = std::max<float>(m, fabs(energy[i*output.height()+j]));
    }
  }
  /*//Grayscale out
  for (unsigned int i=0; i<output.width(); i++) {
    for (unsigned int j=0; j<output.height(); j++) {
      output(i, j, 0) = 255.0*pow(energy[i*output.height()+j]/m, 1/3.0);
      //output(i, j, 1) = 255.0*pow(energy[i*output.height()+j]/m, 1/3.0);
      //output(i, j, 2) = 255.0*pow(energy[i*output.height()+j]/m, 1/3.0);
    }
  }*/

  //Color out
  for (unsigned int i=0; i<output.width(); i++) {
    for (unsigned int j=0; j<output.height(); j++) {
      output(i, j, 0) = image[i*output.height()+j][0];
      output(i, j, 1) = image[i*output.height()+j][1];
      output(i, j, 2) = image[i*output.height()+j][2];
    }
  }
  
  CImg<double> rgb = output.LabtoRGB();
  if (strstr(argv[2], "png"))
    rgb.save_png(argv[2]);
  else if (strstr(argv[2], "jpg"))
    rgb.save_jpeg(argv[2]);
  
  delete [] image;
  return 0;
}
