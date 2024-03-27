/* main.c: Sample of using CXFORM in C
**
** 2003/09/12  Ryan Boller:  Initial version
** 2003/11/11  Ryan Boller:  Added date2es function
** 2004/11/25  Ryan Boller:  Moved date2es and other utility functions to cxform-manual shared lib
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cxform.h"



int main()
{
	int retVal, year, month, day, hour, minute, second;
  int esi;
  long es;
	double jd;
	Vec v_in, v_out0, v_out1, v_out2, v_out3;
	static char inSys[] = "J2000";
	static char outSys0[] = "GEO";
	static char outSys1[] = "GSE";
	static char outSys2[] = "GSM";
	static char outSys3[] = "SM";

  if (0) {
  	year = 1910;
  	month = 1;
  	day = 1;
  	hour = 0;
  	minute = 0;
  	second = 0;

    printf("----------\n");
  	jd = gregorian_calendar_to_jd(year, month, day, hour, minute, second);
  	esi = date2es(year, month, day, hour, minute, second);
    es = date2es(year, month, day, hour, minute, second);
  	
  	printf("Time: %.4d/%.2d/%.2d %.2d:%.2d:%.2d  (JD %f, ES (long): %ld, ES (int): %d)\n\n", year, month, day, hour, minute, second, jd, es, esi);

    //double fracYearIndex = (es+3155803200.0)/157788000.0;
    double fracYearIndex = ((double)es+3155673600.0)/(5.0*365.2*86400.0);
    double fracYear = fmod(fracYearIndex, 1.0);

    printf("fracYear = %f\n", fracYear);
    printf("fracYearIndex = %f\n", fracYearIndex);

    year = 2000;
    month = 1;
    day = 1;
    hour = 0;
    minute = 0;
    second = 0;

    printf("\n----------\n");
    jd = gregorian_calendar_to_jd(year, month, day, hour, minute, second);
    es = date2es(year, month, day, hour, minute, second);
    
    printf("Time: %.4d/%.2d/%.2d %.2d:%.2d:%.2d  (JD %f, ES (long): %ld)\n\n", year, month, day, hour, minute, second, jd, es);
    //fracYearIndex = (es+3155803200.0)/157788000.0;
    fracYearIndex = ((double)es+3155673600.0)/(5.0*365.2*86400.0);
    fracYear = fmod(fracYearIndex, 1.0);

    printf("fracYear = %f\n", fracYear);
    printf("fracYearIndex = %f\n", fracYearIndex);

    year = 1905;
    month = 1;
    day = 1;
    hour = 0;
    minute = 0;
    second = 0;

    printf("\n----------\n");
    jd = gregorian_calendar_to_jd(year, month, day, hour, minute, second);
    es = date2es(year, month, day, hour, minute, second);
    
    printf("Time: %.4d/%.2d/%.2d %.2d:%.2d:%.2d  (JD %f, ES (long): %ld)\n\n", year, month, day, hour, minute, second, jd, es);

    //fracYearIndex = (es+3155803200.0)/157788000.0;
    fracYearIndex = ((double)es+3155673600.0)/(5.0*365.2*86400.0);
    fracYear = fmod(fracYearIndex, 1.0);

    printf("fracYear = %f\n", fracYear);
    printf("fracYearIndex = %f\n", fracYearIndex);

    year = 1900;
    month = 1;
    day = 1;
    hour = 0;
    minute = 0;
    second = 0;

    printf("\n----------\n");
    jd = gregorian_calendar_to_jd(year, month, day, hour, minute, second);
    es = date2es(year, month, day, hour, minute, second);
    
    printf("Time: %.4d/%.2d/%.2d %.2d:%.2d:%.2d  (JD %f, ES (long): %ld)\n\n", year, month, day, hour, minute, second, jd, es);

    //fracYearIndex = (es+3155803200.0)/157788000.0;
    fracYearIndex = ((double)es+3155673600.0)/(5.0*365.2*86400.0);
    fracYear = fmod(fracYearIndex, 1.0);

    printf("fracYear = %f\n", fracYear);
    printf("fracYearIndex = %f\n", fracYearIndex);


    for (year=1900;year < 2030;year = year+5) {
      printf("\n----------\n");
      jd = gregorian_calendar_to_jd(year, month, day, hour, minute, second);
      es = date2es(year, month, day, hour, minute, second);
      
      printf("Time: %.4d/%.2d/%.2d %.2d:%.2d:%.2d  (JD %f, ES (long): %ld)\n\n", year, month, day, hour, minute, second, jd, es);

      //fracYearIndex = (es+3155803200.0)/157788000.0;
      fracYearIndex = ((double)es+3155673600.0)/(5.0*365.2*86400.0);
      fracYear = fmod(fracYearIndex, 1.0);

      printf("fracYear = %f\n", fracYear);
      printf("fracYearIndex = %f\n", fracYearIndex);
    }	
    //return 0;
  }

  year = 2000;
  month = 1;
  day = 1;
  hour = 0;
  minute = 0;
  second = 0;
  jd = gregorian_calendar_to_jd(year, month, day, hour, minute, second);
  es = date2es(year, month, day, hour, minute, second);

	v_in[0] = 1.0;
	v_in[1] = 0.0; 
	v_in[2] = 0.0;

	retVal = cxform(outSys0, outSys2, (double) es, v_in, v_out0);
	//retVal = cxform(inSys, outSys1, es, v_in, v_out1);
	//retVal = cxform(inSys, outSys2, es, v_in, v_out2);
	//retVal = cxform(inSys, outSys3, es, v_in, v_out3);

	if (retVal == 0)
	{
		printf("Input Vector (%s): \t%f\t%f\t%f\n",
			outSys0, v_in[0], v_in[1], v_in[2]);
			
		printf("Output Vector (%s):\t%f\t%f\t%f\n",
			outSys2, v_out0[0], v_out0[1], v_out0[2]);
		//printf("Output Vector (%s):\t%f\t%f\t%f\n",
		//	outSys1, v_out1[0], v_out1[1], v_out1[2]);
		//printf("Output Vector (%s):\t%f\t%f\t%f\n",
		//	outSys2, v_out2[0], v_out2[1], v_out2[2]);
		//printf("Output Vector (%s):\t%f\t%f\t%f\n",
		//	outSys3, v_out3[0], v_out3[1], v_out3[2]);						
	}
	else
		printf("Error during call to cxform\n");
		
	fflush(stdout);
	
	return 0; 
}
