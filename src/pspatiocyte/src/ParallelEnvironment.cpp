#include "ParallelEnvironment.hpp"

//=========================================================================
//   partition Nitem items into Nsect sections
//=========================================================================
void ParallelEnvironment::partition( const int &Nitem, const int &Nsect,
                                     int bgnitem[], int enditem[],
                                     const int &rank, const bool &verbose )
{
   const int itps = Nitem/Nsect;    // least items per section
   const int frac = Nitem%Nsect;

   for( int i=0; i<Nsect; ++i )
   {
      if( i<frac )
      {
         bgnitem[i] =  i   *(itps+1);
         enditem[i] = (i+1)*(itps+1) - 1;
      } else {
         bgnitem[i] =  i   *itps + frac;
         enditem[i] = (i+1)*itps + frac - 1;
      } 
   }

   // check if system length is even
   bool error = false;
   for( int i=0; i<Nsect; ++i )
   {
       const int Npoint = enditem[i] - bgnitem[i] + 1;
       if( Npoint%2 != 0 ) error = true;  // if lengh of axis is odd
   }
   cart.Allreduce( MPI::IN_PLACE, &error, 1, MPI::BOOL, MPI::LOR );
   
   if(error) {
     if (!rank) {
       std::cout << "ERROR: at least one of the specified lattice dimensions" <<
         " results in an odd number after dividing the voxels to" <<
         " subprocesses. This will cause undefined behaviour. Exiting." <<
         std::endl;  
     }
     fout << "ERROR: at least one of the specified lattice dimensions" <<
       " results in an odd number after dividing the voxels to" <<
       " subprocesses. This will cause undefined behaviour. Exiting." <<
       std::endl;  
     MPI::Finalize();
     exit(0);
   }

   if( rank!=0 ) return;

   if( true )
   {
      fout << "  Npoint  = " << setw(4) << Nitem << endl;
      fout << "  Nsector = " << setw(4) << Nsect << endl;
      for( int i=0; i<Nsect; ++i )
      {
         fout << setw(8) << i << "s = ";
         fout << setw(4) << bgnitem[i] + 1 << " to "
              << setw(4) << enditem[i] + 1 << " ("
              << setw(4) << enditem[i] - bgnitem[i] + 1 << "pts)" << endl;
      }
   }
}

