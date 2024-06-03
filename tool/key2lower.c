#include <string.h>
#include <ctype.h>

void key2lower(char *key){
  unsigned int ii;
  for (ii = 0; ii < strlen(key); ii++) key[ii] = tolower(key[ii]);
}
