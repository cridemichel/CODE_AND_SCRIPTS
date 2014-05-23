#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>

#define MC_STEPS 100000
#define BOX_SIZE_X 100
#define BOX_SIZE_Y 100
#define BOX_SIZE_Z 100
#define RADIUS 9.5
struct vector {
		float x;
		float y;
		float z;
	};
float P_x[5000000], P_y[5000000], P_z[5000000];

struct vector P[5000000];
int found_one=0;
////********  main starts here  ********/////
int main( int argc, char *argv[] ){

  	FILE *traj;
	FILE *e2e;
	FILE *ra;
  	int opt;
        float cc;
 	//compute average
	float x,y,z;
	float l,m;
	float comx, comy, comz;
	char a[10000], c[10000], d[10000], e[10000];
	int b,res;

		struct vector PA1_1, PA2_1, PB2_1, PB1_1;
	struct vector bar1_1, bar2_1;

	float dist_1;
	float dist_ist, dist_aver=0.;
	int frame=0;

	FILE *buffer_read;
		
	char input[100]="traj1_fit.pdb", output[100]="e2e.dat";

	FILE *buffer;
	
	char string[100];
        while ((opt = getopt (argc, argv, "i:e:o:")) != -1){

		switch (opt){
			case 'i':
        	        	printf ("Input file: \"%s\"\n", optarg);
				traj=fopen( optarg, "r" );
        	        	break;
			case 'e':
        	        	printf ("Output file: \"%s\"\n", optarg);
				e2e=fopen( optarg, "w" );
        		        break;
			case 'o':
        	        	printf ("Output file: \"%s\"\n", optarg);
				ra=fopen( optarg, "w" );
        	        break;
		//	case 'h':
			//	printf("./nameprogram -i input_trajectori.pdb -e output_end2end -o output_angle");
    		}
  	}




        buffer=fopen("buffer.pdb", "w");
 

	//buffering P positions into a file
	while ( fgets(string, 100, traj) != NULL ){
		if( string[13] == 'P'  ){
			fprintf(buffer, "%s", string);
		}
	}

	fclose(traj);
	fclose(buffer);



        buffer_read=fopen("buffer.pdb", "r");

	while ( !feof(buffer_read) ){
	       	fscanf(buffer_read, "%22c %d %f %f %f %f %f\n", a, &res, &x, &y, &z, &l, &m);

#if 0
	  fscanf(buffer_read, "%d %s %s %d ",  &b, &*c, &*d, &res);
	  fscanf(buffer_read, "%f %f %f ", &x, &y, &z);
	  fscanf(buffer_read, "%f %f\n", &l, &m);
#endif
	   //printf("%f %f %f \n", x, y, z);
	   //printf("%d ", res);

	  if ((res-1)%24==2-1)       { PA1_1.x=x; PA1_1.y=y; PA1_1.z=z; }
	  else if ((res-1)%24==12-1) { PA2_1.x=x; PA2_1.y=y; PA2_1.z=z; }
	  else if ((res-1)%24==14-1) { PB2_1.x=x; PB2_1.y=y; PB2_1.z=z; }
	  else if ((res-1)%24==24-1) { PB1_1.x=x; PB1_1.y=y; PB1_1.z=z; }
	  

	  if((res-1)%24==24-1){
	    bar1_1.x=(PA1_1.x+PB1_1.x)/2.; bar1_1.y=(PA1_1.y+PB1_1.y)/2.; bar1_1.z=(PA1_1.z+PB1_1.z)/2.;
	    bar2_1.x=(PA2_1.x+PB2_1.x)/2.; bar2_1.y=(PA2_1.y+PB2_1.y)/2.; bar2_1.z=(PA2_1.z+PB2_1.z)/2.;
	    
	    // printf("%f %f %f     ", bar1_1.x,bar1_1.y,bar1_1.z);
	    // printf("%f %f %f     \n", bar2_1.x,bar2_1.y,bar2_1.z);
	    
	    dist_ist= pow ( (bar1_1.x-bar2_1.x)*(bar1_1.x-bar2_1.x) + (bar1_1.y-bar2_1.y)*(bar1_1.y-bar2_1.y) + (bar1_1.z-bar2_1.z)*(bar1_1.z-bar2_1.z) , 0.5);
	    
	    frame++;
	    dist_aver = dist_aver+dist_ist;
	    
	    fprintf(e2e, "%d %f %f\n", frame, dist_ist, dist_aver/frame);
	  }
	}
       	printf("\nafter %d steps the average e2e is:  %f\n", frame, dist_aver/frame);

	fclose(e2e);
	//end first section


	//check if e2e of baricentral of P-first group is the same of the averave value
	rewind(buffer_read);
	frame=0;

#if 0
	float dist_final;

	PA1_1.x=0.; PA1_1.y=0.; PA1_1.z=0.;
	PA2_1.x=0.; PA2_1.y=0.; PA2_1.z=0.;
	PB2_1.x=0.; PB2_1.y=0.; PB2_1.z=0.;
	PB1_1.x=0.; PB1_1.y=0.; PB1_1.z=0.;
#if 1
 while ( !feof(buffer_read)){
	       	fscanf(buffer_read, "%22c %d %f %f %f %f %f\n", a, &res, &x, &y, &z, &l, &m);

  		if(res==2)      { PA1_1.x+=x; PA1_1.y+=y; PA1_1.z+=z; }
       		else if(res==12){ PA2_1.x+=x; PA2_1.y+=y; PA2_1.z+=z; }
       	        else if(res==14){ PB2_1.x+=x; PB2_1.y+=y; PB2_1.z+=z; }
       		else if(res==24){ PB1_1.x+=x; PB1_1.y+=y; PB1_1.z+=z; }

     		if(res==24)  frame++;

		//printf("%f %f %f \n", P_x[P_count], P_y[P_count], P_z[P_count]);
	}
#else
	while ( fscanf(buffer_read, "%s", &*a) == 1 ){
		fscanf(buffer_read, "%d %s %s %d ",  &b, &*c, &*d, &res);
	       	fscanf(buffer_read, "%f %f %f ", &x, &y, &z);
	       	fscanf(buffer_read, "%f %f\n", &l, &m);

       		// printf("%f %f %f \n", x, y, z);
       		// printf("%d ", res);

       		if(res==2)      { PA1_1.x+=x; PA1_1.y+=y; PA1_1.z+=z; }
       		else if(res==12){ PA2_1.x+=x; PA2_1.y+=y; PA2_1.z+=z; }
       	        else if(res==14){ PB2_1.x+=x; PB2_1.y+=y; PB2_1.z+=z; }
       		else if(res==24){ PB1_1.x+=x; PB1_1.y+=y; PB1_1.z+=z; }

     		if(res==24)  frame++;
	}
#endif
	PA1_1.x=PA1_1.x/frame; PA1_1.y=PA1_1.y/frame; PA1_1.z=PA1_1.z/frame;
	PA2_1.x=PA2_1.x/frame; PA2_1.y=PA2_1.y/frame; PA2_1.z=PA2_1.z/frame;
	PB2_1.x=PB2_1.x/frame; PB2_1.y=PB2_1.y/frame; PB2_1.z=PB2_1.z/frame;
	PB1_1.x=PB1_1.x/frame; PB1_1.y=PB1_1.y/frame; PB1_1.z=PB1_1.z/frame;

	bar1_1.x=(PA1_1.x+PB1_1.x)/2.; bar1_1.y=(PA1_1.y+PB1_1.y)/2.; bar1_1.z=(PA1_1.z+PB1_1.z)/2.;
	bar2_1.x=(PA2_1.x+PB2_1.x)/2.; bar2_1.y=(PA2_1.y+PB2_1.y)/2.; bar2_1.z=(PA2_1.z+PB2_1.z)/2.;

        dist_final = pow ( (bar1_1.x-bar2_1.x)*(bar1_1.x-bar2_1.x) + (bar1_1.y-bar2_1.y)*(bar1_1.y-bar2_1.y) + (bar1_1.z-bar2_1.z)*(bar1_1.z-bar2_1.z) , 0.5);
  
       	printf("after %d steps e2e of average P is: %f\n\n", frame, dist_final);

	printf("bar1: %f %f %f \n", bar1_1.x,bar1_1.y,bar1_1.z);
	printf("bar2: %f %f %f \n", bar2_1.x,bar2_1.y,bar2_1.z);
	//**** END second section ****//
#endif

	//**** START third section ****//
	//store every P atoms to compute montecarlo
	rewind(buffer_read);
	int P_count=0, k;
#if 1
       while ( !feof(buffer_read)){
	       	fscanf(buffer_read, "%22c %d %f %f %f %f %f\n", a, &res, &x, &y, &z, &l, &m);

		P_x[P_count] = x;
		P_y[P_count] = y;
		P_z[P_count] = z;

		P[P_count].x = x;
		P[P_count].y = y;
		P[P_count].z = z;
		P_count++;
                if (P_count%10000==0) printf("P_count=%d\n", P_count);
		//printf("%f %f %f \n", P_x[P_count], P_y[P_count], P_z[P_count]);
	}

#else
	
	while ( fscanf(buffer_read, "%s", &*a) == 1 ){
	 	fscanf(buffer_read, "%d %s %s %d ",  &b, &*c, &*d, &res);
	       	fscanf(buffer_read, "%f %f %f ", &x, &y, &z);
	       	fscanf(buffer_read, "%f %f\n", &l, &m);

		P_x[P_count] = x;
		P_y[P_count] = y;
		P_z[P_count] = z;

		P[P_count].x = x;
		P[P_count].y = y;
		P[P_count].z = z;
		P_count++;

		//printf("%f %f %f \n", P_x[P_count], P_y[P_count], P_z[P_count]);
	}
#endif	
	//for (k=0; k<P_count; k++)		printf("%f %f %f \n", P_x[k], P_y[k], P_z[k]);

	printf("number of total P atoms stored: %d\n\n", P_count);
	fclose(buffer_read);


	//montecarlo alghoritm begin here
	int i,j, count;
	struct vector bar1, bar2;
	struct vector random[MC_STEPS], d1, d2;
	struct vector v1, v2, u1, u2;
	struct vector w1, w2, q1, q2;
	struct vector base1_A, base1_B, base2_A, base2_B;
	
	char rand[100]="check.dat";
	float phi, phi_new, phi_new_min, l1, l2, l1_new,l2_new,l1_new_min,l2_new_min;
	float rad1, rad2;
	float dist1, p1_x, p1_y;
	float dist2, p2_x, p2_y;
	float rmsd=0., rmsd_min=100000000.;
	float sphere_check=0., ellips_check=0.;

	float l1_aver=0., l2_aver=0., angle_aver=0.;


	for(j=0; j<MC_STEPS; j++){
			random[j].x = (BOX_SIZE_X * drand48() );
			random[j].y = (BOX_SIZE_Y * drand48() );
			random[j].z = (BOX_SIZE_Z * drand48() );
	}
        printf("qui P_count=%d\n", P_count);
        cc=0;
        srand48(145);
	for(i=0; i<P_count; i=i+22){
		rmsd_min=100000000.;
		ellips_check=0.;
		//printf("%f %f %f \n", P[i+10].x,P[i+10].y,P[i+10].z);

		//center of mass of P for each step
		bar1.x = (P[i].x+P[i+21].x)*0.5;
		bar1.y = (P[i].y+P[i+21].y)*0.5;
		bar1.z = (P[i].z+P[i+21].z)*0.5;

		bar2.x = (P[i+10].x+P[i+11].x)*0.5;
		bar2.y = (P[i+10].y+P[i+11].y)*0.5;
		bar2.z = (P[i+10].z+P[i+11].z)*0.5;

		//printf("%f %f %f \n", bar1.x,bar1.y,bar1.z);
		found_one=0;
                comx=comy=comz=0.0;	
		for(k=0; k<22; k++)
		{
                       comx += P[i+k].x;
                       comy += P[i+k].y;
		       comz += P[i+k].z;

                }
                comx /=22.;
                comy /=22.;
                comz /=22.;
		for(j=0; j<MC_STEPS; j++){
			count=0;
			random[j].x = comx+(BOX_SIZE_X * (drand48()-0.5));
			random[j].y = comy+(BOX_SIZE_Y * (drand48()-0.5));
			random[j].z = comz+(BOX_SIZE_Z * (drand48()-0.5));
		
			//printf ( "%f %f %f \n", random.x, random.y, random.z);	

			//compute interesting stuff	
			v1.x = bar1.x - random[j].x;	
			v1.y = bar1.y - random[j].y;
			v1.z = bar1.z - random[j].z;	
			l1 = pow( v1.x*v1.x + v1.y*v1.y + v1.z*v1.z, 0.5 );
			u1.x=v1.x/l1;
			u1.y=v1.y/l1;
			u1.z=v1.z/l1;

			v2.x = bar2.x - random[j].x;
			v2.y = bar2.y - random[j].y;
			v2.z = bar2.z - random[j].z;
			l2 = pow( v2.x*v2.x + v2.y*v2.y + v2.z*v2.z, 0.5 );
			u2.x=v2.x/l2;
			u2.y=v2.y/l2;
			u2.z=v2.z/l2;
	
			phi = (360./6.28319) * acos( (v1.x*v2.x + v1.y*v2.y + v1.z*v2.z)/(l1*l2) );
	
			ellips_check = l1 + l2;
			//printf ("%d %f %f %f \n", i, l1, l2, phi);

			
			//float cylindrical fit
			if( (ellips_check < 36.) && (ellips_check > 34.) && (l1/l2 <1.3) && (l1/l2 >0.7) ){

				//if(phi>=120.)
				//fprintf (ra, "%d %f %f %f \n", i, l1, l2, phi);	
	
			        found_one=1;	
				//find the centers of each basis
				base1_A.x = random[j].x + u1.x*(l1+1.5);
				base1_A.y = random[j].y + u1.y*(l1+1.5);
				base1_A.z = random[j].z + u1.z*(l1+1.5);
		
				base1_B.x = bar1.x - u1.x*(l1+1.5);
				base1_B.y = bar1.y - u1.y*(l1+1.5);
				base1_B.z = bar1.z - u1.z*(l1+1.5);
		
				w1.x = base1_A.x - base1_B.x;
				w1.y = base1_A.y - base1_B.y;
				w1.z = base1_A.z - base1_B.z;
				l1_new = pow( w1.x*w1.x + w1.y*w1.y + w1.z*w1.z , 0.5 );
				q1.x = w1.x/l1_new;
				q1.y = w1.y/l1_new;
				q1.z = w1.z/l1_new;
		
		
				base2_A.x = random[j].x + u2.x*(l2+1.5);
				base2_A.y = random[j].y + u2.y*(l2+1.5);
				base2_A.z = random[j].z + u2.z*(l2+1.5);
			
				base2_B.x = bar2.x - u2.x*(l2+1.5);
				base2_B.y = bar2.y - u2.y*(l2+1.5);
				base2_B.z = bar2.z - u2.z*(l2+1.5);
		
				w2.x = base2_A.x - base2_B.x;
				w2.y = base2_A.y - base2_B.y;
				w2.z = base2_A.z - base2_B.z;
				l2_new = pow( w2.x*w2.x + w2.y*w2.y + w2.z*w2.z , 0.5 );
				q2.x = w2.x/l2_new;
				q2.y = w2.y/l2_new;
				q2.z = w2.z/l2_new;
		
				phi_new = (360./6.28319) * acos( (w1.x*w2.x + w1.y*w2.y + w1.z*w2.z)/(l1_new*l2_new) );	
		
	
			       	//fprintf (ra, "%d %f %f %f %f %f \n", i, l1, l2, l1_new, l2_new, phi);
				//fprintf (ra, "%d %f %f %f\n", i, l1-l1_new, l2-l2_new, phi);
		
				//runtime to compute if point P are inside cylinders
				rmsd = 0.;

				for(k=0; k<22; k++){
					//cylinder 1
					d1.x = P[k+i].x - base1_B.x;
					d1.y = P[k+i].y - base1_B.y;
					d1.z = P[k+i].z - base1_B.z;

					dist1 = pow( d1.x*d1.x + d1.y*d1.y + d1.z*d1.z , 0.5 );
					p1_x = ( d1.x*q1.x + d1.y*q1.y + d1.z*q1.z);
					p1_y = pow( dist1*dist1 - p1_x*p1_x, 0.5 );
		
				 	//cylinder 2
					d2.x = P[k+i].x - base2_B.x;
					d2.y = P[k+i].y - base2_B.y;
					d2.z = P[k+i].z - base2_B.z;
		
					dist2 = pow( d2.x*d2.x + d2.y*d2.y + d2.z*d2.z , 0.5 );
					p2_x = ( d2.x*q2.x + d2.y*q2.y + d2.z*q2.z);
					p2_y = pow( dist2*dist2 - p2_x*p2_x, 0.5 );
		
					rad1 = (p1_y-RADIUS)*(p1_y-RADIUS);
					rad2 = (p2_y-RADIUS)*(p2_y-RADIUS);
		
					if ( rad1 < rad2 ) rmsd = rmsd + rad1;
				 	if ( rad1 >= rad2 ) rmsd = rmsd + rad2;
					//printf("*");
				      }
		
				//fprintf(ra, "%f %f %f %lf\n", l1_new, l2_new, phi_new, rmsd);

				if(rmsd < rmsd_min){
					rmsd_min = rmsd;
					l1_new_min = l1_new;
					l2_new_min = l2_new;
					phi_new_min = phi_new; 
				}

			}

		}
                if (found_one)
                {
		  angle_aver = angle_aver + phi_new_min;
		  l1_aver = l1_aver + l1_new_min;
		  l2_aver = l2_aver + l2_new_min;
                  cc++;
                }
             if ((i/22)%200==0 && cc > 0.0) 
		printf("#found: %f -- %f %f %f %f\n", cc,l1_aver/cc, l2_aver/cc, angle_aver/cc, rmsd_min);
		fprintf(ra, "%d %f %f %f 	%lf %f %f 	%f\n", cc, l1_aver/cc, l2_aver/cc, angle_aver/cc, rmsd_min);
      }
  fclose(ra);
  return 0;
}
