

#include "a3r.h"
#include "vect3d.h"
#include "util.h"
#include <stdio.h>
#include "timelog.h"

//-------------------------------------------------------------------------------------------

A3R_HEADER::A3R_HEADER() 
: count(0)
, r(0)
, count_1(0)
{ 
	strcpy(file_type, "a3r"); 
	strcpy(version, "a"); 
	data_start = sizeof(A3R_HEADER);
}

//-------------------------------------------------------------------------------------------

const int nbuf = 40;
const int perm_header[40] = 
{
  0,   1,  2,  3,
  7,   6,  5,  4,
  11, 10,  9,  8,
  12, 13, 14, 15, 16, 17, 18, 19, 20, 21,
  22, 23, 
  31, 30, 29, 28, 27, 26, 25, 24, 
  35, 34, 33, 32, 
  36, 37, 38, 39
};

//-------------------------------------------------------------------------------------------

int Save_A3R(const char* file_name, Vect3D* start, float* colors, int n, double r) 
{
  FILE *f;
  if (0 == (f = fopen(file_name, "wb")) )     return 0;

	A3R_HEADER header;
	header.count = n;
	header.r = r;

  strcpy(header.version, "d"); 

  if (!is_inverse_byte_order()) // IBM-процессор
  {
    Log.PushMessage("Save_A3R :: Using IBM byte order");
    char header_buf[nbuf];
    char *p = reinterpret_cast<char*>(&header);
    for (int k = 0; k < nbuf; k++)
    {
      header_buf[perm_header[k]] = p[k];
    }
    fwrite(&header_buf, nbuf, 1, f);

    {
	    for (Vect3D* i = start; i != start + n; i++) 
      {
  		  Vect3D::SFLOAT sr = *i;
        char *p = reinterpret_cast<char*>(&sr);
        {
          char q[4] = {p[3], p[2], p[1], p[0]};
          fwrite(&q, 1, 4, f);
        }
        {
          char q[4] = {p[7], p[6], p[5], p[4]};
          fwrite(&q, 1, 4, f);
        }
        {
          char q[4] = {p[11], p[10], p[9], p[8]};
          fwrite(&q, 1, 4, f);
        }
      }
    }
    {
	    for (float* i = colors; i != colors + n; i++) 
      {
  		  float c = *i;
        char *p = reinterpret_cast<char*>(&c);
        {
          char q[4] = {p[3], p[2], p[1], p[0]};
          fwrite(&q, 1, 4, f);
        }
  //      Log.PushMessage("c = ", c);
      }
    }

  }
  else
  {  //Intel
    Log.PushMessage("Save_A3R :: Using inverse byte order");

    fwrite(&header, sizeof(header), 1, f);
    {
	    for (Vect3D* i = start; i != start + n; i++) 
      {
  		  Vect3D::SFLOAT sr = *i; 
        fwrite(&sr,  1, sizeof(Vect3D::SFLOAT), f);
      }
    }
    {
      for (float* i = colors; i != colors + n; i++) 
      {
        fwrite(i, 1, sizeof(float), f);
    //    Log.PushMessage("c = ", *i);

      }
    }
  }
  fclose(f);		

return 1;}

//-------------------------------------------------------------------------------------------

Vect3D* Load_A3R(const char* file_name, int& n, double& r)	
{
  FILE *f;
  if (0 == (f = fopen(file_name, "rb")) )     return 0;

	A3R_HEADER header;

  fread(&header, sizeof(A3R_HEADER),1 , f);
	if(strcmp(header.file_type, "a3r")) return NULL;
	n = header.count;
	r = header.r;

	Vect3D::SFLOAT* buf = new Vect3D::SFLOAT[n];
	Vect3D::SFLOAT* j = buf;

  fread(buf, sizeof(Vect3D::SFLOAT), n, f);
  fclose(f);

	
	Vect3D* start = new Vect3D[n];
	Vect3D* stop = start + n;
	for (Vect3D* i = start; i != stop; i++) *i = *j++;
	delete [] buf;
	
	return start;
}

//-------------------------------------------------------------------------------------------


