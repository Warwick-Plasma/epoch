#include "CommonCommands.h"
#include <stdio.h>
void GetCompositeName(char* Name, char* Class, char*&Composite)
{
    sprintf(Composite,"%s/%s",Class,Name);
}
