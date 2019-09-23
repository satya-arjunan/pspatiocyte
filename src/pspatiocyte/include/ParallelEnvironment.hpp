#ifndef __PARAENV_HPP
#define __PARAENV_HPP

#include <mpi.h>
#include <iomanip>
#include <boost/filesystem.hpp>
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


class ParallelEnvironment
{
public:
    ParallelEnvironment(int argc, char* argv[], const int &Nx,
                        const int &Ny, const int &Nz,
                        const std::string dirname)
    {
        Nx_ = Nx;
        Ny_ = Ny;
        Nz_ = Nz;
        comm = MPI::COMM_WORLD;   // initially allocated processes
        verbose = false;           // diagnostics

        // start MPI
        MPI::Init( argc, argv );

        // 3D cartesian decomposition
        size = comm.Get_size();
        ndims = 3;                // 3 dimensional decomposition
        dims[2] =
        dims[1] =
        dims[0] = 0;              // corrupt without initialization
        MPI::Compute_dims( size, ndims, dims );

        // 3D cartesian communicator
        periods[2] =
        periods[1] =
        periods[0] = false;       // no periodicity
        reorder = true;           // allow processor reordering
        cart = comm.Create_cart( ndims, dims, periods, reorder );
        size = cart.Get_size();
        rank = cart.Get_rank();

        // 3D process coordinates
        cart.Get_coords( rank, ndims, coords );

        // identify neighbor processes
        cart.Shift( 2, -1, inbnd[5], outbnd[5] );  // i<<
        cart.Shift( 2,  1, inbnd[4], outbnd[4] );  // i>>
        cart.Shift( 1, -1, inbnd[3], outbnd[3] );  // j<<
        cart.Shift( 1,  1, inbnd[2], outbnd[2] );  // j>>
        cart.Shift( 0, -1, inbnd[1], outbnd[1] );  // k<<
        cart.Shift( 0,  1, inbnd[0], outbnd[0] );  // k>>

        // prepare parallel files
        char fname[80], fname2[80], fname3[80];
        boost::filesystem::create_directory(dirname);

//        sprintf( fname, "xout%06d", rank );   // bug fix for RICC
        sprintf( fname, "./%s/%06d", dirname.c_str(), rank );
        fout.open( fname );

        sprintf( fname2, "./%s/coordinates_%06d", dirname.c_str(), rank );
        fout2.open( fname2 );

        sprintf( fname3, "./%s/timecourses_%06d", dirname.c_str(), rank );
        fout3.open( fname3 );


        // diagnostics
        if( verbose )
        {
           if( rank==0 )
           {
              int  charlen, vers, subvers;
              char name[80];
              MPI::Get_processor_name( name, charlen );
              MPI::Get_version( vers, subvers );
              fout << "hostname = " << name << endl
                   << "MPI-" << vers << "." << subvers << endl
                   << "number of processes = " << size << endl;
              fout << "decomposition = (" << dims[2] << ","
                                          << dims[1] << ","
                                          << dims[0] << ")" << endl << endl;
           }

           fout << "  rank = " << rank
                << ": cart = (" << coords[2] << ","
                                << coords[1] << ","
                                << coords[0] << ") :" << endl
                << "  i+ = " << inbnd[4] << ">" << rank << ">" << outbnd[4]
                << "  j+ = " << inbnd[2] << ">" << rank << ">" << outbnd[2]
                << "  k+ = " << inbnd[0] << ">" << rank << ">" << outbnd[0]
                << endl
                << "  i- = " << outbnd[5] << "<" << rank << "<" << inbnd[5]
                << "  j- = " << outbnd[3] << "<" << rank << "<" << inbnd[3]
                << "  k- = " << outbnd[1] << "<" << rank << "<" << inbnd[1]
                << endl
                << endl;   // c.f. MPI::PROC_NULL = -2 (Linux)

           fout.flush();
           cart.Barrier();

           if( !topologycheck() ) { cart.Abort( -1 ); }
        }

        // lattice decomposition
        ibgn = new int[dims[2]];   iend = new int[dims[2]];
        jbgn = new int[dims[1]];   jend = new int[dims[1]];
        kbgn = new int[dims[0]];   kend = new int[dims[0]];

        partition( Nx, dims[2], ibgn, iend, rank, verbose );
        partition( Ny, dims[1], jbgn, jend, rank, verbose );
        partition( Nz, dims[0], kbgn, kend, rank, verbose );

        ib = ibgn[coords[2]];   ie = iend[coords[2]];   in = ie - ib + 1;
        jb = jbgn[coords[1]];   je = jend[coords[1]];   jn = je - jb + 1;
        kb = kbgn[coords[0]];   ke = kend[coords[0]];   kn = ke - kb + 1;

        // ib,ie etc. represents global coordinates without halo
        // locally, 0 to in-1 etc. are mapped onto 0+halo to in-1+halo etc.
    }

    ~ParallelEnvironment()
    {
        // terminate MPI
        MPI::Finalize();

        fout.close();
        fout2.close();
        fout3.close();

        delete ibgn, iend;
        delete jbgn, jend;
        delete kbgn, kend;
    }

    unsigned getNx() {
      return Nx_;
    }

    unsigned getNy() {
      return Ny_;
    }

    unsigned getNz() {
      return Nz_;
    }

    int getndims()
    {
        return ndims;
    }

    int getsize()
    {
        return size;
    }

    int getrank()
    {
        return rank;
    }

    int getdim0()
    {
        return dims[0];
    }

    int getdim1()
    {
        return dims[1];
    }

    int getdim2()
    {
        return dims[2];
    }

    int getcoord0()
    {
        return coords[0];
    }

    int getcoord1()
    {
        return coords[1];
    }

    int getcoord2()
    {
        return coords[2];
    }

    int getneighbor_i_plus()
    {
        return outbnd[0];
    }

    int getneighbor_i_minus()
    {
        return inbnd[0];
    }

    int getneighbor_j_plus()
    {
        return outbnd[2];
    }

    int getneighbor_j_minus()
    {
        return inbnd[2];
    }

    int getneighbor_k_plus()
    {
        return outbnd[4];
    }

    int getneighbor_k_minus()
    {
        return inbnd[4];
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
        allcount = new counter[size];
        int           repeat[2] = { 3, 18 };
        MPI::Datatype  types[2] = { MPI::INT, MPI::DOUBLE };
        MPI::Aint     displs[2] = { 0, 3*sizeof(int) };
        MPI::Datatype stype;
        stype = MPI::Datatype::Create_struct( 2, repeat, displs, types );
        stype.Commit();
        cart.Gather( &localcount, 1, stype, allcount, 1, stype, 0 );

       if(rank == 0) {
         std::cout << localcount.atota << std::endl;
       }
        // print statistics
       if(rank==0 && verbose)
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
           for( int i=0; i<size; i++ )
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
           for( int i=0; i<size; i++ )
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
           ave[0] /= size;
           ave[1] /= size;
           ave[2] /= size;
           ave[3] /= size;
           ave[4] /= size;
           ave[5] /= size;
           ave[6] /= size;

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

    int     ndims;   // dimension of cartesian decomposition
    int      size;   // size of cartesian communicator
    int      rank;   // rank of cartesian communicator
    int   dims[3];   // expansion in each axis
    int coords[3];   // process coordinates
    int  inbnd[6];   // inbound processes
    int outbnd[6];   // outbound processes
    unsigned Nx_;
    unsigned Ny_;
    unsigned Nz_;

    bool reorder;    // processor reordering
    bool periods[3]; // periodicity
    bool verbose = false;    // diagnostics

    int *ibgn, *iend;
    int *jbgn, *jend;
    int *kbgn, *kend;

    int ib, ie, in,
        jb, je, jn,
        kb, ke, kn;

    struct counter *allcount;  // performance counter of all processes
    struct counter localcount; // performance counter of local

    void partition( const int &Nitem, const int &Nsect,
                    int bgnitem[], int enditem[],
                    const int &rank, const bool &verbose );

    bool topologycheck( void );
};

#endif /* __PARAENV_HPP */
