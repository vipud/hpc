/* Jake Wise
 * 2018-09-22
 * Program that determines if a year, based on the Gregorian calendar, is a leap year by using algorithm that
 * first divides the year by 4. If it is divible by 4, then divide by 100. 
 * If it is divisble by 100, then check with  * 400. If the number is evenly divisble by 400, 
 * then it is a leap year.
 * Any year is a leap year if it is divisble by 4, but not 100. If it is divisble by 100,
 * it must also be divisble by * 400. 
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

int main() {
  char line[128];
  long yearToTest = 1;
  
  printf("Enter a year to test. Enter '0' to quit.\n");
  while(fgets(line, sizeof line, stdin) != NULL)
    { 
      char *eptr = NULL;
      yearToTest = strtol(line, &eptr, 10);
      if(eptr == NULL || !isspace(*eptr))
	{
	  printf("Invalid input: %s", line);
	  yearToTest = 1;
	  continue;
	}
      if(!yearToTest)
	break;
      
      if(yearToTest%4) {
	printf("Year is not a leap year.\n");
      }
      else if(yearToTest%100) {
	printf("Year is a leap year.\n");
      }
      else if(yearToTest%400) {
	printf("Year is not a leap year.\n");
      }
      else {
	printf("Year is a leap year.\n");
      }
      printf("Enter a year to test. Enter '0' to quit.\n");
    }

  return 0;
}
