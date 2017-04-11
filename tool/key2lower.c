#include <string.h>
#include <ctype.h>

void key2lower(char *key){
  int ii;
  for (ii = 0; ii < strlen(key); ii++) key[ii] = tolower(key[ii]);
}
