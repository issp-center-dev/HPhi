/* HPhi  -  Quantum Lattice Model Simulator */
/* Copyright (C) 2015 The University of Tokyo */

/* This program is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>. */
/*-------------------------------------------------------------
 *[ver.2008.9.1]
 *  多次元配列を確保するマクロ (malloc, free関数 と 引数つき#defineマクロ使用)
 *  使用する場合: int mfint[7] を宣言
 *-------------------------------------------------------------
 * Copyright (C) 2007-2009 Daisuke Tahara. All rights reserved.
 *-------------------------------------------------------------*/


/*=================================================================================================*/
#ifndef HPHI_MFMEMORY_H
#define HPHI_MFMEMORY_H

/*complex type*/
#define c_malloc1(X, N1) \
        X = (double complex*)malloc((N1)*sizeof(double complex));

#define c_malloc2(X, N1, N2) \
        X = (double complex**)malloc((N1)*sizeof(double complex*));\
        for(mfint[0]=0;mfint[0]<(N1);mfint[0]++){\
                X[mfint[0]] = (double complex*)malloc((N2)*sizeof(double complex));\
        }

#define c_malloc3(X, N1, N2, N3) \
        X = (double complex***)malloc((N1)*sizeof(double complex**));\
        for(mfint[0]=0;mfint[0]<(N1);mfint[0]++){\
                X[mfint[0]] = (double complex**)malloc((N2)*sizeof(double complex*));\
                for(mfint[1]=0;mfint[1]<(N2);mfint[1]++){\
                        X[mfint[0]][mfint[1]] = (double complex*)malloc((N3)*sizeof(double complex));\
                }\
        }

#define c_free1(X, N1) \
        free(X);

#define c_free2(X, N1, N2) \
        for(mfint[0]=0;mfint[0]<(N1);mfint[0]++){\
                free(X[mfint[0]]);\
        }\
        free(X);

#define c_free3(X, N1, N2, N3) \
        for(mfint[0]=0;mfint[0]<(N1);mfint[0]++){\
                for(mfint[1]=0;mfint[1]<(N2);mfint[1]++){\
                        free(X[mfint[0]][mfint[1]]);\
                }\
                free(X[mfint[0]]);\
        }\
        free(X);




/*char type*/
#define char_malloc2(X, N1, N2) \
	X = (char**)malloc((N1)*sizeof(char*));\
	for(mfint[0]=0;mfint[0]<(N1);mfint[0]++){\
		X[mfint[0]] = (char*)malloc((N2)*sizeof(char));\
	}

/*double type*/
#define d_malloc1(X, N1) \
	X = (double*)malloc((N1)*sizeof(double));

#define d_malloc2(X, N1, N2) \
	X = (double**)malloc((N1)*sizeof(double*));\
	for(mfint[0]=0;mfint[0]<(N1);mfint[0]++){\
		X[mfint[0]] = (double*)malloc((N2)*sizeof(double));\
	}

#define d_malloc3(X, N1, N2, N3) \
	X = (double***)malloc((N1)*sizeof(double**));\
	for(mfint[0]=0;mfint[0]<(N1);mfint[0]++){\
		X[mfint[0]] = (double**)malloc((N2)*sizeof(double*));\
		for(mfint[1]=0;mfint[1]<(N2);mfint[1]++){\
			X[mfint[0]][mfint[1]] = (double*)malloc((N3)*sizeof(double));\
		}\
	}

#define d_malloc4(X, N1, N2, N3, N4) \
	X = (double****)malloc((N1)*sizeof(double***));\
	for(mfint[0]=0;mfint[0]<(N1);mfint[0]++){\
		X[mfint[0]] = (double***)malloc((N2)*sizeof(double**));\
		for(mfint[1]=0;mfint[1]<(N2);mfint[1]++){\
			X[mfint[0]][mfint[1]] = (double**)malloc((N3)*sizeof(double*));\
			for(mfint[2]=0;mfint[2]<(N3);mfint[2]++){\
				X[mfint[0]][mfint[1]][mfint[2]] = (double*)malloc((N4)*sizeof(double));\
			}\
		}\
	}

#define d_malloc5(X, N1, N2, N3, N4, N5) \
	X = (double*****)malloc((N1)*sizeof(double****));\
	for(mfint[0]=0;mfint[0]<(N1);mfint[0]++){\
		X[mfint[0]] = (double****)malloc((N2)*sizeof(double***));\
		for(mfint[1]=0;mfint[1]<(N2);mfint[1]++){\
			X[mfint[0]][mfint[1]] = (double***)malloc((N3)*sizeof(double**));\
			for(mfint[2]=0;mfint[2]<(N3);mfint[2]++){\
				X[mfint[0]][mfint[1]][mfint[2]] = (double**)malloc((N4)*sizeof(double*));\
				for(mfint[3]=0;mfint[3]<(N4);mfint[3]++){\
					X[mfint[0]][mfint[1]][mfint[2]][mfint[3]] = (double*)malloc((N5)*sizeof(double));\
				}\
			}\
		}\
	}

#define d_malloc6(X, N1, N2, N3, N4, N5, N6) \
	X = (double******)malloc((N1)*sizeof(double*****));\
	for(mfint[0]=0;mfint[0]<(N1);mfint[0]++){\
		X[mfint[0]] = (double*****)malloc((N2)*sizeof(double****));\
		for(mfint[1]=0;mfint[1]<(N2);mfint[1]++){\
			X[mfint[0]][mfint[1]] = (double****)malloc((N3)*sizeof(double***));\
			for(mfint[2]=0;mfint[2]<(N3);mfint[2]++){\
				X[mfint[0]][mfint[1]][mfint[2]] = (double***)malloc((N4)*sizeof(double**));\
				for(mfint[3]=0;mfint[3]<(N4);mfint[3]++){\
					X[mfint[0]][mfint[1]][mfint[2]][mfint[3]] = (double**)malloc((N5)*sizeof(double*));\
					for(mfint[4]=0;mfint[4]<(N5);mfint[4]++){\
						X[mfint[0]][mfint[1]][mfint[2]][mfint[3]][mfint[4]] = (double*)malloc((N6)*sizeof(double));\
					}\
				}\
			}\
		}\
	}

#define d_malloc7(X, N1, N2, N3, N4, N5, N6, N7) \
	X = (double*******)malloc((N1)*sizeof(double******));\
	for(mfint[0]=0;mfint[0]<(N1);mfint[0]++){\
		X[mfint[0]] = (double******)malloc((N2)*sizeof(double*****));\
		for(mfint[1]=0;mfint[1]<(N2);mfint[1]++){\
			X[mfint[0]][mfint[1]] = (double*****)malloc((N3)*sizeof(double****));\
			for(mfint[2]=0;mfint[2]<(N3);mfint[2]++){\
				X[mfint[0]][mfint[1]][mfint[2]] = (double****)malloc((N4)*sizeof(double***));\
				for(mfint[3]=0;mfint[3]<(N4);mfint[3]++){\
					X[mfint[0]][mfint[1]][mfint[2]][mfint[3]] = (double***)malloc((N5)*sizeof(double**));\
					for(mfint[4]=0;mfint[4]<(N5);mfint[4]++){\
						X[mfint[0]][mfint[1]][mfint[2]][mfint[3]][mfint[4]] = (double**)malloc((N6)*sizeof(double*));\
						for(mfint[5]=0;mfint[5]<(N6);mfint[5]++){\
							X[mfint[0]][mfint[1]][mfint[2]][mfint[3]][mfint[4]][mfint[5]] \
							= (double*)malloc((N7)*sizeof(double));\
						}\
					}\
				}\
			}\
		}\
	}

#define d_free1(X, N1) \
	free(X);

#define d_free2(X, N1, N2) \
	for(mfint[0]=0;mfint[0]<(N1);mfint[0]++){\
		free(X[mfint[0]]);\
	}\
	free(X);

#define d_free3(X, N1, N2, N3) \
	for(mfint[0]=0;mfint[0]<(N1);mfint[0]++){\
		for(mfint[1]=0;mfint[1]<(N2);mfint[1]++){\
			free(X[mfint[0]][mfint[1]]);\
		}\
		free(X[mfint[0]]);\
	}\
	free(X);

#define d_free4(X, N1, N2, N3, N4) \
	for(mfint[0]=0;mfint[0]<(N1);mfint[0]++){\
		for(mfint[1]=0;mfint[1]<(N2);mfint[1]++){\
			for(mfint[2]=0;mfint[2]<(N3);mfint[2]++){\
				free(X[mfint[0]][mfint[1]][mfint[2]]);\
			}\
			free(X[mfint[0]][mfint[1]]);\
		}\
		free(X[mfint[0]]);\
	}\
	free(X);

#define d_free5(X, N1, N2, N3, N4, N5) \
	for(mfint[0]=0;mfint[0]<(N1);mfint[0]++){\
		for(mfint[1]=0;mfint[1]<(N2);mfint[1]++){\
			for(mfint[2]=0;mfint[2]<(N3);mfint[2]++){\
				for(mfint[3]=0;mfint[3]<(N4);mfint[3]++){\
					free(X[mfint[0]][mfint[1]][mfint[2]][mfint[3]]);\
				}\
				free(X[mfint[0]][mfint[1]][mfint[2]]);\
			}\
			free(X[mfint[0]][mfint[1]]);\
		}\
		free(X[mfint[0]]);\
	}\
	free(X);

#define d_free6(X, N1, N2, N3, N4, N5, N6) \
	for(mfint[0]=0;mfint[0]<(N1);mfint[0]++){\
		for(mfint[1]=0;mfint[1]<(N2);mfint[1]++){\
			for(mfint[2]=0;mfint[2]<(N3);mfint[2]++){\
				for(mfint[3]=0;mfint[3]<(N4);mfint[3]++){\
					for(mfint[4]=0;mfint[4]<(N5);mfint[4]++){\
						free(X[mfint[0]][mfint[1]][mfint[2]][mfint[3]][mfint[4]]);\
					}\
					free(X[mfint[0]][mfint[1]][mfint[2]][mfint[3]]);\
				}\
				free(X[mfint[0]][mfint[1]][mfint[2]]);\
			}\
			free(X[mfint[0]][mfint[1]]);\
		}\
		free(X[mfint[0]]);\
	}\
	free(X);

#define d_free7(X, N1, N2, N3, N4, N5, N6, N7) \
	for(mfint[0]=0;mfint[0]<(N1);mfint[0]++){\
		for(mfint[1]=0;mfint[1]<(N2);mfint[1]++){\
			for(mfint[2]=0;mfint[2]<(N3);mfint[2]++){\
				for(mfint[3]=0;mfint[3]<(N4);mfint[3]++){\
					for(mfint[4]=0;mfint[4]<(N5);mfint[4]++){\
						for(mfint[5]=0;mfint[5]<(N6);mfint[5]++){\
							free(X[mfint[0]][mfint[1]][mfint[2]][mfint[3]][mfint[4]][mfint[5]]);\
						}\
						free(X[mfint[0]][mfint[1]][mfint[2]][mfint[3]][mfint[4]]);\
					}\
					free(X[mfint[0]][mfint[1]][mfint[2]][mfint[3]]);\
				}\
				free(X[mfint[0]][mfint[1]][mfint[2]]);\
			}\
			free(X[mfint[0]][mfint[1]]);\
		}\
		free(X[mfint[0]]);\
	}\
	free(X);

/*unsigned int*/
#define ui_malloc1(X, N1) \
	X = (unsigned int *)malloc((N1)*sizeof(unsigned int));
/*long int*/
#define li_malloc1(X, N1) \
	X = (long int *)malloc((N1)*sizeof(long int));

#define li_malloc2(X, N1, N2) \
	X = (long int**)malloc((N1)*sizeof(long int*));\
	for(mfint[0]=0;mfint[0]<(N1);mfint[0]++){\
		X[mfint[0]] = (long int*)malloc((N2)*sizeof(long int));\
	}



/*long unsigned int*/
#define lui_malloc1(X, N1) \
	X = (long unsigned int *)malloc((N1)*sizeof(long unsigned int));

/*int type*/
#define i_malloc1(X, N1) \
	X = (int*)malloc((N1)*sizeof(int));

#define i_malloc2(X, N1, N2) \
	X = (int**)malloc((N1)*sizeof(int*));\
	for(mfint[0]=0;mfint[0]<(N1);mfint[0]++){\
		X[mfint[0]] = (int*)malloc((N2)*sizeof(int));\
	}

#define li_malloc2(X, N1, N2) \
	X = (long int**)malloc((N1)*sizeof(long int*));\
	for(mfint[0]=0;mfint[0]<(N1);mfint[0]++){\
		X[mfint[0]] = (long int*)malloc((N2)*sizeof(long int));\
	}

#define i_malloc3(X, N1, N2, N3) \
	X = (int***)malloc((N1)*sizeof(int**));\
	for(mfint[0]=0;mfint[0]<(N1);mfint[0]++){\
		X[mfint[0]] = (int**)malloc((N2)*sizeof(int*));\
		for(mfint[1]=0;mfint[1]<(N2);mfint[1]++){\
			X[mfint[0]][mfint[1]] = (int*)malloc((N3)*sizeof(int));\
		}\
	}

#define i_malloc4(X, N1, N2, N3, N4) \
	X = (int****)malloc((N1)*sizeof(int***));\
	for(mfint[0]=0;mfint[0]<(N1);mfint[0]++){\
		X[mfint[0]] = (int***)malloc((N2)*sizeof(int**));\
		for(mfint[1]=0;mfint[1]<(N2);mfint[1]++){\
			X[mfint[0]][mfint[1]] = (int**)malloc((N3)*sizeof(int*));\
			for(mfint[2]=0;mfint[2]<(N3);mfint[2]++){\
				X[mfint[0]][mfint[1]][mfint[2]] = (int*)malloc((N4)*sizeof(int));\
			}\
		}\
	}

#define i_malloc5(X, N1, N2, N3, N4, N5) \
	X = (int*****)malloc((N1)*sizeof(int****));\
	for(mfint[0]=0;mfint[0]<(N1);mfint[0]++){\
		X[mfint[0]] = (int****)malloc((N2)*sizeof(int***));\
		for(mfint[1]=0;mfint[1]<(N2);mfint[1]++){\
			X[mfint[0]][mfint[1]] = (int***)malloc((N3)*sizeof(int**));\
			for(mfint[2]=0;mfint[2]<(N3);mfint[2]++){\
				X[mfint[0]][mfint[1]][mfint[2]] = (int**)malloc((N4)*sizeof(int*));\
				for(mfint[3]=0;mfint[3]<(N4);mfint[3]++){\
					X[mfint[0]][mfint[1]][mfint[2]][mfint[3]] = (int*)malloc((N5)*sizeof(int));\
				}\
			}\
		}\
	}

#define i_malloc6(X, N1, N2, N3, N4, N5, N6) \
	X = (int******)malloc((N1)*sizeof(int*****));\
	for(mfint[0]=0;mfint[0]<(N1);mfint[0]++){\
		X[mfint[0]] = (int*****)malloc((N2)*sizeof(int****));\
		for(mfint[1]=0;mfint[1]<(N2);mfint[1]++){\
			X[mfint[0]][mfint[1]] = (int****)malloc((N3)*sizeof(int***));\
			for(mfint[2]=0;mfint[2]<(N3);mfint[2]++){\
				X[mfint[0]][mfint[1]][mfint[2]] = (int***)malloc((N4)*sizeof(int**));\
				for(mfint[3]=0;mfint[3]<(N4);mfint[3]++){\
					X[mfint[0]][mfint[1]][mfint[2]][mfint[3]] = (int**)malloc((N5)*sizeof(int*));\
					for(mfint[4]=0;mfint[4]<(N5);mfint[4]++){\
						X[mfint[0]][mfint[1]][mfint[2]][mfint[3]][mfint[4]] = (int*)malloc((N6)*sizeof(int));\
					}\
				}\
			}\
		}\
	}

/*unsigned int*/
#define ui_free(X, N1) \
	free(X);
#define i_free1(X, N1) \
	free(X);

#define i_free2(X, N1, N2) \
	for(mfint[0]=0;mfint[0]<(N1);mfint[0]++){\
		free(X[mfint[0]]);\
	}\
	free(X);

#define i_free3(X, N1, N2, N3) \
	for(mfint[0]=0;mfint[0]<(N1);mfint[0]++){\
		for(mfint[1]=0;mfint[1]<(N2);mfint[1]++){\
			free(X[mfint[0]][mfint[1]]);\
		}\
		free(X[mfint[0]]);\
	}\
	free(X);

#define i_free4(X, N1, N2, N3, N4) \
	for(mfint[0]=0;mfint[0]<(N1);mfint[0]++){\
		for(mfint[1]=0;mfint[1]<(N2);mfint[1]++){\
			for(mfint[2]=0;mfint[2]<(N3);mfint[2]++){\
				free(X[mfint[0]][mfint[1]][mfint[2]]);\
			}\
			free(X[mfint[0]][mfint[1]]);\
		}\
		free(X[mfint[0]]);\
	}\
	free(X);

#define i_free5(X, N1, N2, N3, N4, N5) \
	for(mfint[0]=0;mfint[0]<(N1);mfint[0]++){\
		for(mfint[1]=0;mfint[1]<(N2);mfint[1]++){\
			for(mfint[2]=0;mfint[2]<(N3);mfint[2]++){\
				for(mfint[3]=0;mfint[3]<(N4);mfint[3]++){\
					free(X[mfint[0]][mfint[1]][mfint[2]][mfint[3]]);\
				}\
				free(X[mfint[0]][mfint[1]][mfint[2]]);\
			}\
			free(X[mfint[0]][mfint[1]]);\
		}\
		free(X[mfint[0]]);\
	}\
	free(X);

#define i_free6(X, N1, N2, N3, N4, N5, N6) \
	for(mfint[0]=0;mfint[0]<(N1);mfint[0]++){\
		for(mfint[1]=0;mfint[1]<(N2);mfint[1]++){\
			for(mfint[2]=0;mfint[2]<(N3);mfint[2]++){\
				for(mfint[3]=0;mfint[3]<(N4);mfint[3]++){\
					for(mfint[4]=0;mfint[4]<(N5);mfint[4]++){\
						free(X[mfint[0]][mfint[1]][mfint[2]][mfint[3]][mfint[4]]);\
					}\
					free(X[mfint[0]][mfint[1]][mfint[2]][mfint[3]]);\
				}\
				free(X[mfint[0]][mfint[1]][mfint[2]]);\
			}\
			free(X[mfint[0]][mfint[1]]);\
		}\
		free(X[mfint[0]]);\
	}\
	free(X);

#endif /* HPHI_MFMEMORY_H */
