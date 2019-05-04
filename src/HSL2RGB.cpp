/*
    A spartial constraint tumour growth model
    Copyright (C) 2018 Timon Heide (timon.heide@icr.ac.uk)
                     & Kate Chkhaidze (Kate.Chkhaidze@icr.ac.uk)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "HSL2RGB.hpp"
//

unsigned char GetRValue(unsigned color) {
	return (unsigned char)((color>>16)&0xFF);
}

unsigned char GetGValue(unsigned color) {
	return (unsigned char)((color>>8)&0xFF);
}

unsigned char GetBValue(unsigned color) {
	return (unsigned char)(color&0xFF);
}

static inline unsigned RGB(unsigned char r,unsigned char g,unsigned char b){
	unsigned color = ((unsigned)r<<16) | ((unsigned)g<<8) | b;
	return color;
}

// This is a subfunction of HSLtoRGB
static void HSLtoRGB_Subfunction(unsigned& c, const float& temp1,
																 const float& temp2, const float& temp3)
{
	if((temp3 * 6) < 1)
		c = (unsigned)((temp2 + (temp1 - temp2)*6*temp3)*100);
	else
		if((temp3 * 2) < 1)
			c = (unsigned)(temp1*100);
		else
			if((temp3 * 3) < 2)
				c = (unsigned)((temp2 + (temp1 - temp2)*(.66666 - temp3)*6)*100);
			else
				c = (unsigned)(temp2*100);
	return;
}


// This function extracts the hue, saturation, and luminance from "color"
// and places these values in h, s, and l respectively.
void RGBtoHSL(unsigned color, unsigned& h, unsigned& s, unsigned& l) {
	unsigned r = (unsigned)GetRValue(color);
	unsigned g = (unsigned)GetGValue(color);
	unsigned b = (unsigned)GetBValue(color);

	float r_percent = ((float)r) / 255;
	float g_percent = ((float)g) / 255;
	float b_percent = ((float)b) / 255;

	float max_color = 0;
	if((r_percent >= g_percent) && (r_percent >= b_percent))
	{
		max_color = r_percent;
	}
	if((g_percent >= r_percent) && (g_percent >= b_percent))
		max_color = g_percent;
	if((b_percent >= r_percent) && (b_percent >= g_percent))
		max_color = b_percent;

	float min_color = 0;
	if((r_percent <= g_percent) && (r_percent <= b_percent))
		min_color = r_percent;
	if((g_percent <= r_percent) && (g_percent <= b_percent))
		min_color = g_percent;
	if((b_percent <= r_percent) && (b_percent <= g_percent))
		min_color = b_percent;

	float L = 0;
	float S = 0;
	float H = 0;

	L = (max_color + min_color) / 2;

	if(max_color == min_color) {
		S = 0;
		H = 0;
	} else {
		if(L < .50) {
			S = (max_color - min_color) / (max_color + min_color);
		} else {
			S = (max_color - min_color) / (2 - max_color - min_color);
		}

		if(max_color == r_percent) {
			H = (g_percent - b_percent)/(max_color - min_color);
		}

		if(max_color == g_percent) {
			H = 2 + (b_percent - r_percent)/(max_color - min_color);
		}

		if(max_color == b_percent) {
			H = 4 + (r_percent - g_percent)/(max_color - min_color);
		}
	}
	s = (unsigned)(S*100);
	l = (unsigned)(L*100);
	H = H*60;
	if(H < 0)
		H += 360;
	h = (unsigned)H;
}

// This function converts the "color" object to the equivalent RGB values of
// the hue, saturation, and luminance passed as h, s, and l respectively
unsigned HSLtoRGB(const unsigned& h, const unsigned& s,
	                    const unsigned& l)
{
	unsigned r = 0;
	unsigned g = 0;
	unsigned b = 0;

	float L = ((float)l) / 100;
	float S = ((float)s) / 100;
	float H = ((float)h) / 360;

	if(s == 0) {
		r = l;
		g = l;
		b = l;
	} else {
		float temp1 = 0;
		if(L < .50) {
			temp1 = L*(1 + S);
		} else {
			temp1 = L + S - (L*S);
		}

		float temp2 = 2*L - temp1;

		float temp3 = 0;
		for(int i = 0 ; i < 3 ; i++)
		{
			switch(i) {
			case 0: // red
				temp3 = H + (1.f / 3.f);
				if(temp3 > 1)
					temp3 -= 1;
				HSLtoRGB_Subfunction(r,temp1,temp2,temp3);
				break;
			case 1: // green
				temp3 = H;
				HSLtoRGB_Subfunction(g,temp1,temp2,temp3);
				break;
			case 2: // blue
				temp3 = H - (1.f / 3.f);
				if(temp3 < 0)
					temp3 += 1;
				HSLtoRGB_Subfunction(b,temp1,temp2,temp3);
				break;
			default:
				{

				}
			}
		}
	}
	r = (unsigned)((((float)r) / 100) * 255);
	g = (unsigned)((((float)g) / 100) * 255);
	b = (unsigned)((((float)b) / 100) * 255);
	return RGB(r, g, b);
}
