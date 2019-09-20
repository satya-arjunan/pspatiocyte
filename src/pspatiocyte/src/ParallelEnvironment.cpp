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
   
   if( error )
   {
       std::cout << "ERROR: odd number lattice will behave illegally. Exiting." << endl;
       fout << "ERROR: odd number lattice will behave illegally. Exiting." << endl;
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

//=========================================================================
//   check process topology of 3D cartesian communicator
//=========================================================================
#define testin( pos, dir ) ( (allin[pos][dir]!=NE) && (allout[allin[pos][dir]][dir]!=pos) )

#define testout( pos, dir )  ( (allout[pos][dir]!=NE) && (allin[allout[pos][dir]][dir]!=pos) )

#define preamble( i ) ( fout << "  CartError: rank = " << i )
#define log_ip( i ) ( preamble( i ) << ": i+ inconsistent" << endl )
#define log_im( i ) ( preamble( i ) << ": i- inconsistent" << endl )
#define log_jp( i ) ( preamble( i ) << ": j+ inconsistent" << endl )
#define log_jm( i ) ( preamble( i ) << ": j- inconsistent" << endl )
#define log_kp( i ) ( preamble( i ) << ": k+ inconsistent" << endl )
#define log_km( i ) ( preamble( i ) << ": k- inconsistent" << endl )

bool ParallelEnvironment::topologycheck( void )
{
//   const int NE = -2;       // neighbor not exist
   const int NE = MPI::PROC_NULL;       // neighbor not exist

   int allcrd[size][3];
   int allin [size][6];
   int allout[size][6];

   cart.Gather( coords, 3, MPI::INT, allcrd, 3, MPI::INT, 0 );
   cart.Gather(  inbnd, 6, MPI::INT, allin,  6, MPI::INT, 0 );
   cart.Gather( outbnd, 6, MPI::INT, allout, 6, MPI::INT, 0 );

//   bool err = false;
   int err = -1;

   if( rank==0 )
   {
      for( int i=0; i<size; ++i )
      {
         // check inbound processes
         if( testin( i, 0 ) )  { log_im( i ); err = 1; }
         if( testin( i, 1 ) )  { log_ip( i ); err = 1; }
         if( testin( i, 2 ) )  { log_jm( i ); err = 1; }
         if( testin( i, 3 ) )  { log_jp( i ); err = 1; }
         if( testin( i, 4 ) )  { log_km( i ); err = 1; }
         if( testin( i, 5 ) )  { log_kp( i ); err = 1; }

         // check outbound processes
         if( testout( i, 0 ) ) { log_ip( i ); err = 1; }
         if( testout( i, 1 ) ) { log_im( i ); err = 1; }
         if( testout( i, 2 ) ) { log_jp( i ); err = 1; }
         if( testout( i, 3 ) ) { log_jm( i ); err = 1; }
         if( testout( i, 4 ) ) { log_kp( i ); err = 1; }
         if( testout( i, 5 ) ) { log_km( i ); err = 1; }
      }
      fout.flush();
   }

//   cart.Bcast( &err, 1, MPI::BOOL, 0 );
   cart.Bcast( &err, 1, MPI::INT, 0 );

//   return err ? false : true;
   return err ? -1 : 1;
}
#undef testin
#undef testout
#undef preamble
#undef log_ip
#undef log_im
#undef log_jp
#undef log_jm
#undef log_kp
#undef log_km

