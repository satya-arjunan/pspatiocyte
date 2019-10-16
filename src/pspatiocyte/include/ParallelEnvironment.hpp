#ifndef __PARAENV_HPP
#define __PARAENV_HPP

#include <mpi.h>
#include <iomanip>
#include <boost/filesystem.hpp>
#include "Vector.hpp"
#include "Common.hpp"
using namespace std;

struct counter
{
    // number of callings
    int ndiff;
    int ninde;
    int ninfl;
    // accumurated time
    double atota;
    double adiff;
    double acalc;
    double acomm;
    double apack;
    double ampic;
    double aupck;
    double ainde;
    double ainfl;
    // lap time
    double ltota;
    double ldiff;
    double lcalc;
    double lcomm;
    double lpack;
    double lmpic;
    double lupck;
    double linde;
    double linfl;
};


class ParallelEnvironment {
public:
  ParallelEnvironment(int argc, char* argv[], const int &GNx, const int &GNy,
                      const int &GNz, const std::string dirname) {
    GNx_ = GNx;
    GNy_ = GNy;
    GNz_ = GNz;
    comm = MPI::COMM_WORLD;   // initially allocated processes
    verbose_ = false;           // diagnostics
    // start MPI
    MPI::Init( argc, argv );
    // 3D cartesian decomposition
    size_ = comm.Get_size();
    ndims_ = 3;                // 3 dimensional decomposition
    dims_[2] =
    dims_[1] =
    dims_[0] = 0;              // corrupt without initialization
    MPI::Compute_dims(size_, ndims_, dims_);
    // 3D cartesian communicator
    periodicity_[2] =
    periodicity_[1] =
    periodicity_[0] = false;       // no periodicity
    reorder_ = true;           // allow processor reorder_ing
    cart = comm.Create_cart( ndims_, dims_, periodicity_, reorder_ );
    size_ = cart.Get_size();
    rank_ = cart.Get_rank();
    // 3D process coordinates
    cart.Get_coords(rank_, ndims_, coords_);
    // identify neighbor processes
    cart.Shift( 2, -1, inbound_[5], outbound_[5] );  // i<<
    cart.Shift( 2,  1, inbound_[4], outbound_[4] );  // i>>
    cart.Shift( 1, -1, inbound_[3], outbound_[3] );  // j<<
    cart.Shift( 1,  1, inbound_[2], outbound_[2] );  // j>>
    cart.Shift( 0, -1, inbound_[1], outbound_[1] );  // k<<
    cart.Shift( 0,  1, inbound_[0], outbound_[0] );  // k>>

    // prepare parallel files
    char fname[80], fname2[80], fname3[80];
    boost::filesystem::create_directory(dirname);

     //sprintf( fname, "xout%06d", rank_ );   // bug fix for RICC
    sprintf( fname, "./%s/%06d", dirname.c_str(), rank_ );
    fout.open( fname );

    sprintf( fname2, "./%s/coordinates_%06d", dirname.c_str(), rank_ );
    fout2.open( fname2 );

    sprintf( fname3, "./%s/timecourses_%06d", dirname.c_str(), rank_ );
    fout3.open( fname3 ); 
    
    if(!rank_) {
      int  charlen, vers, subvers;
      char name[80];
      MPI::Get_processor_name( name, charlen );
      MPI::Get_version( vers, subvers );
      fout << "hostname = " << name << endl
        << "MPI-" << vers << "." << subvers << endl
        << "number of processes = " << size_ << endl;
      fout << "decomposition = (" << dims_[2] << ","
        << dims_[1] << "," << dims_[0] << ")" << endl << endl;
    } 
    fout << "  rank_ = " << rank_ << ": cart = ("
      << coords_[2] << "," << coords_[1] << "," << coords_[0] << ") :" << endl
      << "  i+ = " << inbound_[4] << ">" << rank_ << ">" << outbound_[4]
      << "  j+ = " << inbound_[2] << ">" << rank_ << ">" << outbound_[2]
      << "  k+ = " << inbound_[0] << ">" << rank_ << ">" << outbound_[0]
      << endl
      << "  i- = " << outbound_[5] << "<" << rank_ << "<" << inbound_[5]
      << "  j- = " << outbound_[3] << "<" << rank_ << "<" << inbound_[3]
      << "  k- = " << outbound_[1] << "<" << rank_ << "<" << inbound_[1]
      << endl << endl;
   fout.flush();
   cart.Barrier();

   if(!topologycheck()) {
     cart.Abort( -1 );
   } 
   
   // lattice decomposition
   ibegin_ = new int[dims_[2]];   iend_ = new int[dims_[2]];
   jbegin_ = new int[dims_[1]];   jend_ = new int[dims_[1]];
   kbegin_ = new int[dims_[0]];   kend_ = new int[dims_[0]];
   fout << "dims:" << dims_[2] << " " << dims_[1] << " " << dims_[0] <<
     std::endl;

   partition( GNx, dims_[2], ibegin_, iend_, rank_, verbose_ );
   partition( GNy, dims_[1], jbegin_, jend_, rank_, verbose_ );
   partition( GNz, dims_[0], kbegin_, kend_, rank_, verbose_ );

   ib = ibegin_[coords_[2]];   ie = iend_[coords_[2]];   in = ie - ib + 1;
   jb = jbegin_[coords_[1]];   je = jend_[coords_[1]];   jn = je - jb + 1;
   kb = kbegin_[coords_[0]];   ke = kend_[coords_[0]];   kn = ke - kb + 1;

   fout << "ib:" << ib << " ie:" << ie << " in:" << in << std::endl;
   fout << "jb:" << jb << " je:" << je << " jn:" << jn << std::endl;
   fout << "kb:" << kb << " ke:" << ke << " kn:" << kn << std::endl;
   // ib,ie etc. represents global coordinates without halo
   // locally, 0 to in-1 etc. are mapped onto 0+halo to in-1+halo etc.
  }

  ~ParallelEnvironment() {
    // terminate MPI
    MPI::Finalize();

    fout.close();
    fout2.close();
    fout3.close();

    delete ibegin_, iend_;
    delete jbegin_, jend_;
    delete kbegin_, kend_;
  }

  Vector<unsigned> get_global_dimensions() {
    return Vector<unsigned>(GNx_, GNy_, GNz_);
  }

  Vector<unsigned> get_local_min() {
    return Vector<unsigned>(ib, jb, kb);
  }

  Vector<unsigned> get_local_max() {
    return Vector<unsigned>(ie, je, ke);
  }

  unsigned get_global_seed(unsigned seed) {
    getcart().Bcast(&seed, 1 , MPI_INT, 0);
    return seed;
  }

  int getndims()
  {
      return ndims_;
  }

  int getsize()
  {
      return size_;
  }

  int getrank()
  {
      return rank_;
  }

  int getdim0()
  {
      return dims_[0];
  }

  int getdim1()
  {
      return dims_[1];
  }

  int getdim2()
  {
      return dims_[2];
  }

  int getcoord0()
  {
      return coords_[0];
  }

  int getcoord1()
  {
      return coords_[1];
  }

  int getcoord2()
  {
      return coords_[2];
  }

  int getneighbor_i_plus()
  {
      return outbound_[0];
  }

  int getneighbor_i_minus()
  {
      return inbound_[0];
  }

  int getneighbor_j_plus()
  {
      return outbound_[2];
  }

  int getneighbor_j_minus()
  {
      return inbound_[2];
  }

  int getneighbor_k_plus()
  {
      return outbound_[4];
  }

  int getneighbor_k_minus()
  {
      return inbound_[4];
  }

  int getispan()
  {
      return in;
  }

  int getjspan()
  {
      return jn;
  }

  int getkspan()
  {
      return kn;
  }

  MPI::Cartcomm &getcart()
  {
      return cart;
  }

  void starttimer()
  {
    localcount.atota = 0;
    localcount.ltota = MPI::Wtime();
    /*
      // initialize counters
      localcount.ndiff =
      localcount.ninde =
      localcount.ninfl = 0;
      localcount.atota =
      localcount.adiff =
      localcount.acomm =
      localcount.apack =
      localcount.ampic =
      localcount.aupck =
      localcount.ainde =
      localcount.ainfl = 0.0;
      */
  }
  void starttimerdiff()
  {
    /*
      localcount.ndiff++;
      localcount.ldiff = MPI::Wtime();
      */
  }
  void stoptimerdiff()
  {
      //localcount.adiff += MPI::Wtime() - localcount.ldiff;
  }

  void starttimercomm()
  {
      //localcount.lcomm = MPI::Wtime();
  }
  void stoptimercomm()
  {
      //localcount.acomm += MPI::Wtime() - localcount.lcomm;
  }

  void starttimerpack()
  {
      //localcount.lpack = MPI::Wtime();
  }
  void stoptimerpack()
  {
      //localcount.apack += MPI::Wtime() - localcount.lpack;
  }

  void starttimermpic()
  {
      //localcount.lmpic = MPI::Wtime();
  }
  void stoptimermpic()
  {
      //localcount.ampic += MPI::Wtime() - localcount.lmpic;
  }

  void starttimerupck()
  {
      //localcount.lupck = MPI::Wtime();
  }
  void stoptimerupck()
  {
      //localcount.aupck += MPI::Wtime() - localcount.lupck;
  }

  void starttimerinde()
  {
      /*localcount.ninde++;
      localcount.linde = MPI::Wtime();
      */
  }
  void stoptimerinde()
  {
      //localcount.ainde += MPI::Wtime() - localcount.linde;
  }

  void starttimerinfl()
  {
      //localcount.ninfl++;
      //localcount.linfl = MPI::Wtime();
  }
  void stoptimerinfl()
  {
      //localcount.ainfl += MPI::Wtime() - localcount.linfl;
  }

  void stoptimer()
  {
      // stop all counters
      localcount.atota += MPI::Wtime() - localcount.ltota;

      // print local counters
      localcount.acalc = localcount.adiff - localcount.acomm;
      fout << right;
      fout << endl << "Elapsed time   = "
           << setw(12) << localcount.atota << " sec  ";
      fout << endl << "|- diffusion   = "
           << setw(12) << localcount.adiff << " sec  "
           << setw(12) << localcount.ndiff << " calls";
      fout << endl << "|  |- calc     = "
           << setw(12) << localcount.acalc << " sec  ";
      fout << endl << "|  |- comm     = "
           << setw(12) << localcount.acomm << " sec  ";
      fout << endl << "|     |- pack  = "
           << setw(12) << localcount.apack << " sec  ";
      fout << endl << "|     |- mpi   = "
           << setw(12) << localcount.ampic << " sec  ";
      fout << endl << "|     |- upck  = "
           << setw(12) << localcount.aupck << " sec  ";
      fout << endl << "|- independent = "
           << setw(12) << localcount.ainde << " sec  "
           << setw(12) << localcount.ninde << " calls";
      fout << endl << "|- influenced  = "
           << setw(12) << localcount.ainfl << " sec  "
           << setw(12) << localcount.ninfl << " calls" << endl;

      // gather counters to root
      allcount = new counter[size_];
      int           repeat[2] = { 3, 18 };
      MPI::Datatype  types[2] = { MPI::INT, MPI::DOUBLE };
      MPI::Aint     displs[2] = { 0, 3*sizeof(int) };
      MPI::Datatype stype;
      stype = MPI::Datatype::Create_struct( 2, repeat, displs, types );
      stype.Commit();
      cart.Gather( &localcount, 1, stype, allcount, 1, stype, 0 );

     if(rank_ == 0) {
       std::cout << localcount.atota << std::endl;
     }
      // print statistics
     if(rank_==0 && verbose_)
      {
         fout << endl << "==============================";
         fout << endl << "=====       results      =====";
         fout << endl << "==============================" << endl;
         fout << endl << right << "Statistics:"
              << setw(10) << "elapse"
              << setw(10) << "diffusion"
              << setw(10) << "calc"
              << setw(10) << "comm"
              << setw(10) << "pack"
              << setw(10) << "mpi"
              << setw(10) << "unpack" << endl;
         #ifdef PRINT_NODETIMES
         for( int i=0; i<size_; i++ )
         {
            fout << setw(11) << i << setw(10) << allcount[i].atota
                                  << setw(10) << allcount[i].adiff
                                  << setw(10) << allcount[i].acalc
                                  << setw(10) << allcount[i].acomm
                                  << setw(10) << allcount[i].apack
                                  << setw(10) << allcount[i].ampic
                                  << setw(10) << allcount[i].aupck << endl;
         }
         #endif

         double max[7], min[7], ave[7];
         for( int i=0; i<7; i++ )
         {
            max[i] = 0.0;
            min[i] = 1.0e8;
            ave[i] = 0.0;
         }
         for( int i=0; i<size_; i++ )
         {
           // update maximum value
           if( allcount[i].atota > max[0] ) max[0] = allcount[i].atota;
           if( allcount[i].adiff > max[1] ) max[1] = allcount[i].adiff;
           if( allcount[i].acalc > max[2] ) max[2] = allcount[i].acalc;
           if( allcount[i].acomm > max[3] ) max[3] = allcount[i].acomm;
           if( allcount[i].apack > max[4] ) max[4] = allcount[i].apack;
           if( allcount[i].ampic > max[5] ) max[5] = allcount[i].ampic;
           if( allcount[i].aupck > max[6] ) max[6] = allcount[i].aupck;

           // update minimum value
           if( min[0] > allcount[i].atota ) min[0] = allcount[i].atota;
           if( min[1] > allcount[i].adiff ) min[1] = allcount[i].adiff;
           if( min[2] > allcount[i].acalc ) min[2] = allcount[i].acalc;
           if( min[3] > allcount[i].acomm ) min[3] = allcount[i].acomm;
           if( min[4] > allcount[i].apack ) min[4] = allcount[i].apack;
           if( min[5] > allcount[i].ampic ) min[5] = allcount[i].ampic;
           if( min[6] > allcount[i].aupck ) min[6] = allcount[i].aupck;

           // accumurate value
           ave[0] += allcount[i].atota;
           ave[1] += allcount[i].adiff;
           ave[2] += allcount[i].acalc;
           ave[3] += allcount[i].acomm;
           ave[4] += allcount[i].apack;
           ave[5] += allcount[i].ampic;
           ave[6] += allcount[i].aupck;
         }
         // take average
         ave[0] /= size_;
         ave[1] /= size_;
         ave[2] /= size_;
         ave[3] /= size_;
         ave[4] /= size_;
         ave[5] /= size_;
         ave[6] /= size_;

         fout << endl;
         fout << setw(11) << "max"
              << setw(10) << max[0] << setw(10) << max[1]
              << setw(10) << max[2] << setw(10) << max[3]
              << setw(10) << max[4] << setw(10) << max[5]
              << setw(10) << max[6] << endl;
         fout << setw(11) << "min"
              << setw(10) << min[0] << setw(10) << min[1]
              << setw(10) << min[2] << setw(10) << min[3]
              << setw(10) << min[4] << setw(10) << min[5]
              << setw(10) << min[6] << endl;
         fout << setw(11) << "ave"
              << setw(10) << ave[0] << setw(10) << ave[1]
              << setw(10) << ave[2] << setw(10) << ave[3]
              << setw(10) << ave[4] << setw(10) << ave[5]
              << setw(10) << ave[6] << endl;
/*
         // CSV format
         fout << endl << "##### CSV format #####" << endl
              << setw(15) << "Elapse-ave,"
              << setw(15) << "Elapse-min,"
              << setw(15) << "Elapse-max,"
              << setw(15) << "Diffusion-ave,"
              << setw(15) << "Diffusion-min,"
              << setw(15) << "Diffusion-max,"
              << setw(15) << "Calc-ave,"
              << setw(15) << "Calc-min,"
              << setw(15) << "Calc-max,"
              << setw(15) << "Comm-ave,"
              << setw(15) << "Comm-min,"
              << setw(15) << "Comm-max,"
              << setw(15) << "Pack-ave,"
              << setw(15) << "Pack-min,"
              << setw(15) << "Pack-max,"
              << setw(15) << "MPI-ave,"
              << setw(15) << "MPI-min,"
              << setw(15) << "MPI-max,"
              << setw(15) << "Unpack-ave,"
              << setw(15) << "Unpack-min,"
              << setw(15) << "Unpack-max," << endl;
         fout << setw(15) << ave[0] << ","
              << setw(15) << min[0] << ","
              << setw(15) << max[0] << ","
              << setw(15) << ave[1] << ","
              << setw(15) << min[1] << ","
              << setw(15) << max[1] << ","
              << setw(15) << ave[2] << ","
              << setw(15) << min[2] << ","
              << setw(15) << max[2] << ","
              << setw(15) << ave[3] << ","
              << setw(15) << min[3] << ","
              << setw(15) << max[3] << ","
              << setw(15) << ave[4] << ","
              << setw(15) << min[4] << ","
              << setw(15) << max[4] << ","
              << setw(15) << ave[5] << ","
              << setw(15) << min[5] << ","
              << setw(15) << max[5] << ","
              << setw(15) << ave[6] << ","
              << setw(15) << min[6] << ","
              << setw(15) << max[6] << "," << endl;
*/
      }
  }


private:
    MPI::Intracomm  comm;    // default communicator
    MPI::Cartcomm   cart;    // cartesian communicator

    int     ndims_;   // dimension of cartesian decomposition
    int      size_;   // size_ of cartesian communicator
    int      rank_;   // rank_ of cartesian communicator
    int   dims_[3];   // expansion in each axis
    int coords_[3];   // process coordinates
    int  inbound_[6];   // inbound processes
    int outbound_[6];   // outbound processes
    unsigned GNx_;
    unsigned GNy_;
    unsigned GNz_;

    bool reorder_;    // processor reorder_ing
    bool periodicity_[3]; // periodicity
    bool verbose_ = false;    // diagnostics

    int *ibegin_, *iend_;
    int *jbegin_, *jend_;
    int *kbegin_, *kend_;

    int ib, ie, in,
        jb, je, jn,
        kb, ke, kn;

    struct counter *allcount;  // performance counter of all processes
    struct counter localcount; // performance counter of local

    void partition(const int &Nitem, const int &Nsect,
                    int bgnitem[], int enditem[],
                    const int &rank, const bool &verbose);

    bool topologycheck( void );
};

#endif /* __PARAENV_HPP */
