# 1 "interpolation.c"
# 1 "<built-in>"
# 1 "<command-line>"
# 31 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 32 "<command-line>" 2
# 1 "interpolation.c"
static char help[] = "Interpolation - fdf-curvIB ";

# 1 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscpf.h" 1





# 1 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscvec.h" 1
# 9 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscvec.h"
# 1 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscis.h" 1






# 1 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h" 1
# 14 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
# 1 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscconf.h" 1
# 15 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h" 2
# 1 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscfix.h" 1
# 16 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h" 2
# 88 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
# 1 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscversion.h" 1
# 89 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h" 2
# 108 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
# 1 "/sw/eb/sw/impi/2018.5.288-iccifort-2019.5.281/include64/mpi.h" 1
# 76 "/sw/eb/sw/impi/2018.5.288-iccifort-2019.5.281/include64/mpi.h"
# 1 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/stdint.h" 1 3 4
# 9 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/stdint.h" 3 4
# 1 "/usr/include/stdint.h" 1 3 4
# 25 "/usr/include/stdint.h" 3 4
# 1 "/usr/include/features.h" 1 3 4
# 375 "/usr/include/features.h" 3 4
# 1 "/usr/include/sys/cdefs.h" 1 3 4
# 392 "/usr/include/sys/cdefs.h" 3 4
# 1 "/usr/include/bits/wordsize.h" 1 3 4
# 393 "/usr/include/sys/cdefs.h" 2 3 4
# 376 "/usr/include/features.h" 2 3 4
# 399 "/usr/include/features.h" 3 4
# 1 "/usr/include/gnu/stubs.h" 1 3 4
# 10 "/usr/include/gnu/stubs.h" 3 4
# 1 "/usr/include/gnu/stubs-64.h" 1 3 4
# 11 "/usr/include/gnu/stubs.h" 2 3 4
# 400 "/usr/include/features.h" 2 3 4
# 26 "/usr/include/stdint.h" 2 3 4
# 1 "/usr/include/bits/wchar.h" 1 3 4
# 22 "/usr/include/bits/wchar.h" 3 4
# 1 "/usr/include/bits/wordsize.h" 1 3 4
# 23 "/usr/include/bits/wchar.h" 2 3 4
# 27 "/usr/include/stdint.h" 2 3 4
# 1 "/usr/include/bits/wordsize.h" 1 3 4
# 28 "/usr/include/stdint.h" 2 3 4
# 36 "/usr/include/stdint.h" 3 4

# 36 "/usr/include/stdint.h" 3 4
typedef signed char int8_t;
typedef short int int16_t;
typedef int int32_t;

typedef long int int64_t;







typedef unsigned char uint8_t;
typedef unsigned short int uint16_t;

typedef unsigned int uint32_t;



typedef unsigned long int uint64_t;
# 65 "/usr/include/stdint.h" 3 4
typedef signed char int_least8_t;
typedef short int int_least16_t;
typedef int int_least32_t;

typedef long int int_least64_t;






typedef unsigned char uint_least8_t;
typedef unsigned short int uint_least16_t;
typedef unsigned int uint_least32_t;

typedef unsigned long int uint_least64_t;
# 90 "/usr/include/stdint.h" 3 4
typedef signed char int_fast8_t;

typedef long int int_fast16_t;
typedef long int int_fast32_t;
typedef long int int_fast64_t;
# 103 "/usr/include/stdint.h" 3 4
typedef unsigned char uint_fast8_t;

typedef unsigned long int uint_fast16_t;
typedef unsigned long int uint_fast32_t;
typedef unsigned long int uint_fast64_t;
# 119 "/usr/include/stdint.h" 3 4
typedef long int intptr_t;


typedef unsigned long int uintptr_t;
# 134 "/usr/include/stdint.h" 3 4
typedef long int intmax_t;
typedef unsigned long int uintmax_t;
# 10 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/stdint.h" 2 3 4
# 77 "/sw/eb/sw/impi/2018.5.288-iccifort-2019.5.281/include64/mpi.h" 2
# 170 "/sw/eb/sw/impi/2018.5.288-iccifort-2019.5.281/include64/mpi.h"

# 170 "/sw/eb/sw/impi/2018.5.288-iccifort-2019.5.281/include64/mpi.h"
typedef int MPI_Datatype;
# 354 "/sw/eb/sw/impi/2018.5.288-iccifort-2019.5.281/include64/mpi.h"
typedef int MPI_Comm;




typedef int MPI_Group;



typedef int MPI_Win;







typedef struct ADIOI_FileD *MPI_File;



typedef int MPI_Op;
# 438 "/sw/eb/sw/impi/2018.5.288-iccifort-2019.5.281/include64/mpi.h"
typedef enum MPIR_Win_flavor {
    MPI_WIN_FLAVOR_CREATE = 1,
    MPI_WIN_FLAVOR_ALLOCATE = 2,
    MPI_WIN_FLAVOR_DYNAMIC = 3,
    MPI_WIN_FLAVOR_SHARED = 4
} MPIR_Win_flavor_t;


typedef enum MPIR_Win_model {
    MPI_WIN_SEPARATE = 1,
    MPI_WIN_UNIFIED = 2
} MPIR_Win_model_t;





typedef enum MPIR_Topo_type { MPI_GRAPH=1, MPI_CART=2, MPI_DIST_GRAPH=3 } MPIR_Topo_type;
# 468 "/sw/eb/sw/impi/2018.5.288-iccifort-2019.5.281/include64/mpi.h"
typedef void (MPI_Handler_function) ( MPI_Comm *, int *, ... );
typedef int (MPI_Comm_copy_attr_function)(MPI_Comm, int, void *, void *,
       void *, int *);
typedef int (MPI_Comm_delete_attr_function)(MPI_Comm, int, void *, void *);
typedef int (MPI_Type_copy_attr_function)(MPI_Datatype, int, void *, void *,
       void *, int *);
typedef int (MPI_Type_delete_attr_function)(MPI_Datatype, int, void *, void *);
typedef int (MPI_Win_copy_attr_function)(MPI_Win, int, void *, void *, void *,
      int *);
typedef int (MPI_Win_delete_attr_function)(MPI_Win, int, void *, void *);

typedef void (MPI_Comm_errhandler_function)(MPI_Comm *, int *, ...);
typedef void (MPI_File_errhandler_function)(MPI_File *, int *, ...);
typedef void (MPI_Win_errhandler_function)(MPI_Win *, int *, ...);

typedef MPI_Comm_errhandler_function MPI_Comm_errhandler_fn;
typedef MPI_File_errhandler_function MPI_File_errhandler_fn;
typedef MPI_Win_errhandler_function MPI_Win_errhandler_fn;
# 496 "/sw/eb/sw/impi/2018.5.288-iccifort-2019.5.281/include64/mpi.h"
typedef int MPI_Errhandler;
# 517 "/sw/eb/sw/impi/2018.5.288-iccifort-2019.5.281/include64/mpi.h"
typedef int MPI_Request;


typedef int MPI_Message;


typedef void (MPI_User_function) ( void *, void *, int *, MPI_Datatype * );


typedef int (MPI_Copy_function) ( MPI_Comm, int, void *, void *, void *, int * );
typedef int (MPI_Delete_function) ( MPI_Comm, int, void *, void * );
# 597 "/sw/eb/sw/impi/2018.5.288-iccifort-2019.5.281/include64/mpi.h"
enum MPIR_Combiner_enum {
    MPI_COMBINER_NAMED = 1,
    MPI_COMBINER_DUP = 2,
    MPI_COMBINER_CONTIGUOUS = 3,
    MPI_COMBINER_VECTOR = 4,
    MPI_COMBINER_HVECTOR_INTEGER = 5,
    MPI_COMBINER_HVECTOR = 6,
    MPI_COMBINER_INDEXED = 7,
    MPI_COMBINER_HINDEXED_INTEGER = 8,
    MPI_COMBINER_HINDEXED = 9,
    MPI_COMBINER_INDEXED_BLOCK = 10,
    MPI_COMBINER_STRUCT_INTEGER = 11,
    MPI_COMBINER_STRUCT = 12,
    MPI_COMBINER_SUBARRAY = 13,
    MPI_COMBINER_DARRAY = 14,
    MPI_COMBINER_F90_REAL = 15,
    MPI_COMBINER_F90_COMPLEX = 16,
    MPI_COMBINER_F90_INTEGER = 17,
    MPI_COMBINER_RESIZED = 18,
    MPI_COMBINER_HINDEXED_BLOCK = 19
};


typedef int MPI_Info;
# 647 "/sw/eb/sw/impi/2018.5.288-iccifort-2019.5.281/include64/mpi.h"
typedef long MPI_Aint;
typedef int MPI_Fint;
typedef long long MPI_Count;
# 666 "/sw/eb/sw/impi/2018.5.288-iccifort-2019.5.281/include64/mpi.h"
typedef long long MPI_Offset;







typedef struct MPI_Status {
    int count_lo;
    int count_hi_and_cancelled;
    int MPI_SOURCE;
    int MPI_TAG;
    int MPI_ERROR;
} MPI_Status;


struct MPIR_T_enum_s;
struct MPIR_T_cvar_handle_s;
struct MPIR_T_pvar_handle_s;
struct MPIR_T_pvar_session_s;

typedef struct MPIR_T_enum_s * MPI_T_enum;
typedef struct MPIR_T_cvar_handle_s * MPI_T_cvar_handle;
typedef struct MPIR_T_pvar_handle_s * MPI_T_pvar_handle;
typedef struct MPIR_T_pvar_session_s * MPI_T_pvar_session;


extern struct MPIR_T_pvar_handle_s * const MPI_T_PVAR_ALL_HANDLES;
# 703 "/sw/eb/sw/impi/2018.5.288-iccifort-2019.5.281/include64/mpi.h"
typedef enum MPIR_T_verbosity_t {


    MPIX_T_VERBOSITY_INVALID = 0,


    MPI_T_VERBOSITY_USER_BASIC = 221,
    MPI_T_VERBOSITY_USER_DETAIL,
    MPI_T_VERBOSITY_USER_ALL,

    MPI_T_VERBOSITY_TUNER_BASIC,
    MPI_T_VERBOSITY_TUNER_DETAIL,
    MPI_T_VERBOSITY_TUNER_ALL,

    MPI_T_VERBOSITY_MPIDEV_BASIC,
    MPI_T_VERBOSITY_MPIDEV_DETAIL,
    MPI_T_VERBOSITY_MPIDEV_ALL
} MPIR_T_verbosity_t;

typedef enum MPIR_T_bind_t {


    MPIX_T_BIND_INVALID = 0,


    MPI_T_BIND_NO_OBJECT = 9700,
    MPI_T_BIND_MPI_COMM,
    MPI_T_BIND_MPI_DATATYPE,
    MPI_T_BIND_MPI_ERRHANDLER,
    MPI_T_BIND_MPI_FILE,
    MPI_T_BIND_MPI_GROUP,
    MPI_T_BIND_MPI_OP,
    MPI_T_BIND_MPI_REQUEST,
    MPI_T_BIND_MPI_WIN,
    MPI_T_BIND_MPI_MESSAGE,
    MPI_T_BIND_MPI_INFO
} MPIR_T_bind_t;

typedef enum MPIR_T_scope_t {


    MPIX_T_SCOPE_INVALID = 0,


    MPI_T_SCOPE_CONSTANT = 60438,
    MPI_T_SCOPE_READONLY,
    MPI_T_SCOPE_LOCAL,
    MPI_T_SCOPE_GROUP,
    MPI_T_SCOPE_GROUP_EQ,
    MPI_T_SCOPE_ALL,
    MPI_T_SCOPE_ALL_EQ
} MPIR_T_scope_t;

typedef enum MPIR_T_pvar_class_t {


    MPIX_T_PVAR_CLASS_INVALID = 0,


    MPIR_T_PVAR_CLASS_FIRST = 240,
    MPI_T_PVAR_CLASS_STATE = MPIR_T_PVAR_CLASS_FIRST,
    MPI_T_PVAR_CLASS_LEVEL,
    MPI_T_PVAR_CLASS_SIZE,
    MPI_T_PVAR_CLASS_PERCENTAGE,
    MPI_T_PVAR_CLASS_HIGHWATERMARK,
    MPI_T_PVAR_CLASS_LOWWATERMARK,
    MPI_T_PVAR_CLASS_COUNTER,
    MPI_T_PVAR_CLASS_AGGREGATE,
    MPI_T_PVAR_CLASS_TIMER,
    MPI_T_PVAR_CLASS_GENERIC,
    MPIR_T_PVAR_CLASS_LAST,
    MPIR_T_PVAR_CLASS_NUMBER = MPIR_T_PVAR_CLASS_LAST - MPIR_T_PVAR_CLASS_FIRST
} MPIR_T_pvar_class_t;
# 825 "/sw/eb/sw/impi/2018.5.288-iccifort-2019.5.281/include64/mpi.h"
extern MPI_Fint * MPI_F_STATUS_IGNORE;
extern MPI_Fint * MPI_F_STATUSES_IGNORE;
extern int * const MPI_UNWEIGHTED;
extern int * const MPI_WEIGHTS_EMPTY;
# 842 "/sw/eb/sw/impi/2018.5.288-iccifort-2019.5.281/include64/mpi.h"
typedef struct {
    MPI_Fint count_lo;
    MPI_Fint count_hi_and_cancelled;
    MPI_Fint MPI_SOURCE;
    MPI_Fint MPI_TAG;
    MPI_Fint MPI_ERROR;
} MPI_F08_Status;

extern MPI_F08_Status MPIR_F08_MPI_STATUS_IGNORE_OBJ;
extern MPI_F08_Status MPIR_F08_MPI_STATUSES_IGNORE_OBJ[1];
extern int MPIR_F08_MPI_IN_PLACE;
extern int MPIR_F08_MPI_BOTTOM;


extern MPI_F08_Status *MPI_F08_STATUS_IGNORE;
extern MPI_F08_Status *MPI_F08_STATUSES_IGNORE;
# 866 "/sw/eb/sw/impi/2018.5.288-iccifort-2019.5.281/include64/mpi.h"
typedef int (MPI_Grequest_cancel_function)(void *, int);
typedef int (MPI_Grequest_free_function)(void *);
typedef int (MPI_Grequest_query_function)(void *, MPI_Status *);
typedef int (MPIX_Grequest_poll_function)(void *, MPI_Status *);
typedef int (MPIX_Grequest_wait_function)(int, void **, double, MPI_Status *);
# 996 "/sw/eb/sw/impi/2018.5.288-iccifort-2019.5.281/include64/mpi.h"
typedef int (MPI_Datarep_conversion_function)(void *, MPI_Datatype, int,
             void *, MPI_Offset, void *);
typedef int (MPI_Datarep_extent_function)(MPI_Datatype datatype, MPI_Aint *,
                      void *);
# 1016 "/sw/eb/sw/impi/2018.5.288-iccifort-2019.5.281/include64/mpi.h"
int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
             MPI_Comm comm) ;
int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
             MPI_Comm comm, MPI_Status *status) ;
int MPI_Get_count(const MPI_Status *status, MPI_Datatype datatype, int *count);
int MPI_Bsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm) ;
int MPI_Ssend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm) ;
int MPI_Rsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm) ;
int MPI_Buffer_attach(void *buffer, int size);
int MPI_Buffer_detach(void *buffer_addr, int *size);
int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm, MPI_Request *request) ;
int MPI_Ibsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
               MPI_Comm comm, MPI_Request *request) ;
int MPI_Issend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
               MPI_Comm comm, MPI_Request *request) ;
int MPI_Irsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
               MPI_Comm comm, MPI_Request *request) ;
int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
              MPI_Comm comm, MPI_Request *request) ;
int MPI_Wait(MPI_Request *request, MPI_Status *status);
int MPI_Test(MPI_Request *request, int *flag, MPI_Status *status);
int MPI_Request_free(MPI_Request *request);
int MPI_Waitany(int count, MPI_Request array_of_requests[], int *indx, MPI_Status *status);
int MPI_Testany(int count, MPI_Request array_of_requests[], int *indx, int *flag,
                MPI_Status *status);
int MPI_Waitall(int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]);
int MPI_Testall(int count, MPI_Request array_of_requests[], int *flag,
                MPI_Status array_of_statuses[]);
int MPI_Waitsome(int incount, MPI_Request array_of_requests[], int *outcount,
                 int array_of_indices[], MPI_Status array_of_statuses[]);
int MPI_Testsome(int incount, MPI_Request array_of_requests[], int *outcount,
                 int array_of_indices[], MPI_Status array_of_statuses[]);
int MPI_Iprobe(int source, int tag, MPI_Comm comm, int *flag, MPI_Status *status);
int MPI_Probe(int source, int tag, MPI_Comm comm, MPI_Status *status);
int MPI_Cancel(MPI_Request *request);
int MPI_Test_cancelled(const MPI_Status *status, int *flag);
int MPI_Send_init(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                  MPI_Comm comm, MPI_Request *request) ;
int MPI_Bsend_init(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                   MPI_Comm comm, MPI_Request *request) ;
int MPI_Ssend_init(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                   MPI_Comm comm, MPI_Request *request) ;
int MPI_Rsend_init(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                   MPI_Comm comm, MPI_Request *request) ;
int MPI_Recv_init(void *buf, int count, MPI_Datatype datatype, int source, int tag,
                  MPI_Comm comm, MPI_Request *request) ;
int MPI_Start(MPI_Request *request);
int MPI_Startall(int count, MPI_Request array_of_requests[]);
int MPI_Sendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, int dest,
                 int sendtag, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                 int source, int recvtag, MPI_Comm comm, MPI_Status *status)
                 ;
int MPI_Sendrecv_replace(void *buf, int count, MPI_Datatype datatype, int dest,
                         int sendtag, int source, int recvtag, MPI_Comm comm,
                         MPI_Status *status) ;
int MPI_Type_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype *newtype);
int MPI_Type_vector(int count, int blocklength, int stride, MPI_Datatype oldtype,
                    MPI_Datatype *newtype);
int MPI_Type_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype,
                     MPI_Datatype *newtype);
int MPI_Type_indexed(int count, const int *array_of_blocklengths,
                     const int *array_of_displacements, MPI_Datatype oldtype,
                     MPI_Datatype *newtype);
int MPI_Type_hindexed(int count, const int *array_of_blocklengths,
                      const MPI_Aint *array_of_displacements, MPI_Datatype oldtype,
                      MPI_Datatype *newtype);
int MPI_Type_struct(int count, const int *array_of_blocklengths,
                    const MPI_Aint *array_of_displacements,
                    const MPI_Datatype *array_of_types, MPI_Datatype *newtype);
int MPI_Address(const void *location, MPI_Aint *address);
int MPI_Type_extent(MPI_Datatype datatype, MPI_Aint *extent);
int MPI_Type_size(MPI_Datatype datatype, int *size);
int MPI_Type_lb(MPI_Datatype datatype, MPI_Aint *displacement);
int MPI_Type_ub(MPI_Datatype datatype, MPI_Aint *displacement);
int MPI_Type_commit(MPI_Datatype *datatype);
int MPI_Type_free(MPI_Datatype *datatype);
int MPI_Get_elements(const MPI_Status *status, MPI_Datatype datatype, int *count);
int MPI_Pack(const void *inbuf, int incount, MPI_Datatype datatype, void *outbuf,
             int outsize, int *position, MPI_Comm comm) ;
int MPI_Unpack(const void *inbuf, int insize, int *position, void *outbuf, int outcount,
               MPI_Datatype datatype, MPI_Comm comm) ;
int MPI_Pack_size(int incount, MPI_Datatype datatype, MPI_Comm comm, int *size);
int MPI_Barrier(MPI_Comm comm);
int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)
              ;
int MPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
               int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
               ;
int MPI_Gatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                const int *recvcounts, const int *displs, MPI_Datatype recvtype, int root,
                MPI_Comm comm)
                ;
int MPI_Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
                ;
int MPI_Scatterv(const void *sendbuf, const int *sendcounts, const int *displs,
                 MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                 int root, MPI_Comm comm)
                 ;
int MPI_Allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                  int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                  ;
int MPI_Allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                   const int *recvcounts, const int *displs, MPI_Datatype recvtype, MPI_Comm comm)
                   ;
int MPI_Alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                 int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                 ;
int MPI_Alltoallv(const void *sendbuf, const int *sendcounts, const int *sdispls,
                  MPI_Datatype sendtype, void *recvbuf, const int *recvcounts,
                  const int *rdispls, MPI_Datatype recvtype, MPI_Comm comm)
                  ;
int MPI_Alltoallw(const void *sendbuf, const int sendcounts[], const int sdispls[],
                  const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
                  const int rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm);
int MPI_Exscan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
               MPI_Op op, MPI_Comm comm)
               ;
int MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
               MPI_Op op, int root, MPI_Comm comm)
               ;
int MPI_Op_create(MPI_User_function *user_fn, int commute, MPI_Op *op);
int MPI_Op_free(MPI_Op *op);
int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                  MPI_Op op, MPI_Comm comm)
                  ;
int MPI_Reduce_scatter(const void *sendbuf, void *recvbuf, const int recvcounts[],
                       MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                       ;
int MPI_Scan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
             MPI_Comm comm)
             ;
int MPI_Group_size(MPI_Group group, int *size);
int MPI_Group_rank(MPI_Group group, int *rank);
int MPI_Group_translate_ranks(MPI_Group group1, int n, const int ranks1[], MPI_Group group2,
                              int ranks2[]);
int MPI_Group_compare(MPI_Group group1, MPI_Group group2, int *result);
int MPI_Comm_group(MPI_Comm comm, MPI_Group *group);
int MPI_Group_union(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup);
int MPI_Group_intersection(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup);
int MPI_Group_difference(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup);
int MPI_Group_incl(MPI_Group group, int n, const int ranks[], MPI_Group *newgroup);
int MPI_Group_excl(MPI_Group group, int n, const int ranks[], MPI_Group *newgroup);
int MPI_Group_range_incl(MPI_Group group, int n, int ranges[][3], MPI_Group *newgroup);
int MPI_Group_range_excl(MPI_Group group, int n, int ranges[][3], MPI_Group *newgroup);
int MPI_Group_free(MPI_Group *group);
int MPI_Comm_size(MPI_Comm comm, int *size);
int MPI_Comm_rank(MPI_Comm comm, int *rank);
int MPI_Comm_compare(MPI_Comm comm1, MPI_Comm comm2, int *result);
int MPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm);
int MPI_Comm_dup_with_info(MPI_Comm comm, MPI_Info info, MPI_Comm *newcomm);
int MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm);
int MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm);
int MPI_Comm_free(MPI_Comm *comm);
int MPI_Comm_test_inter(MPI_Comm comm, int *flag);
int MPI_Comm_remote_size(MPI_Comm comm, int *size);
int MPI_Comm_remote_group(MPI_Comm comm, MPI_Group *group);
int MPI_Intercomm_create(MPI_Comm local_comm, int local_leader, MPI_Comm peer_comm,
                         int remote_leader, int tag, MPI_Comm *newintercomm);
int MPI_Intercomm_merge(MPI_Comm intercomm, int high, MPI_Comm *newintracomm);
int MPI_Keyval_create(MPI_Copy_function *copy_fn, MPI_Delete_function *delete_fn,
                      int *keyval, void *extra_state);
int MPI_Keyval_free(int *keyval);
int MPI_Attr_put(MPI_Comm comm, int keyval, void *attribute_val);
int MPI_Attr_get(MPI_Comm comm, int keyval, void *attribute_val, int *flag);
int MPI_Attr_delete(MPI_Comm comm, int keyval);
int MPI_Topo_test(MPI_Comm comm, int *status);
int MPI_Cart_create(MPI_Comm comm_old, int ndims, const int dims[], const int periods[],
                    int reorder, MPI_Comm *comm_cart);
int MPI_Dims_create(int nnodes, int ndims, int dims[]);
int MPI_Graph_create(MPI_Comm comm_old, int nnodes, const int indx[], const int edges[],
                     int reorder, MPI_Comm *comm_graph);
int MPI_Graphdims_get(MPI_Comm comm, int *nnodes, int *nedges);
int MPI_Graph_get(MPI_Comm comm, int maxindex, int maxedges, int indx[], int edges[]);
int MPI_Cartdim_get(MPI_Comm comm, int *ndims);
int MPI_Cart_get(MPI_Comm comm, int maxdims, int dims[], int periods[], int coords[]);
int MPI_Cart_rank(MPI_Comm comm, const int coords[], int *rank);
int MPI_Cart_coords(MPI_Comm comm, int rank, int maxdims, int coords[]);
int MPI_Graph_neighbors_count(MPI_Comm comm, int rank, int *nneighbors);
int MPI_Graph_neighbors(MPI_Comm comm, int rank, int maxneighbors, int neighbors[]);
int MPI_Cart_shift(MPI_Comm comm, int direction, int disp, int *rank_source, int *rank_dest);
int MPI_Cart_sub(MPI_Comm comm, const int remain_dims[], MPI_Comm *newcomm);
int MPI_Cart_map(MPI_Comm comm, int ndims, const int dims[], const int periods[], int *newrank);
int MPI_Graph_map(MPI_Comm comm, int nnodes, const int indx[], const int edges[], int *newrank);
int MPI_Get_processor_name(char *name, int *resultlen);
int MPI_Get_version(int *version, int *subversion);
int MPI_Get_library_version(char *version, int *resultlen);
int MPI_Errhandler_create(MPI_Handler_function *function, MPI_Errhandler *errhandler);
int MPI_Errhandler_set(MPI_Comm comm, MPI_Errhandler errhandler);
int MPI_Errhandler_get(MPI_Comm comm, MPI_Errhandler *errhandler);
int MPI_Errhandler_free(MPI_Errhandler *errhandler);
int MPI_Error_string(int errorcode, char *string, int *resultlen);
int MPI_Error_class(int errorcode, int *errorclass);
double MPI_Wtime(void);
double MPI_Wtick(void);
int MPI_Init(int *argc, char ***argv);
int MPI_Finalize(void);
int MPI_Initialized(int *flag);
int MPI_Abort(MPI_Comm comm, int errorcode);



int MPI_Pcontrol(const int level, ...);
int MPIR_Dup_fn(MPI_Comm oldcomm, int keyval, void *extra_state, void *attribute_val_in,
               void *attribute_val_out, int *flag);


int MPI_Close_port(const char *port_name);
int MPI_Comm_accept(const char *port_name, MPI_Info info, int root, MPI_Comm comm,
                    MPI_Comm *newcomm);
int MPI_Comm_connect(const char *port_name, MPI_Info info, int root, MPI_Comm comm,
                     MPI_Comm *newcomm);
int MPI_Comm_disconnect(MPI_Comm *comm);
int MPI_Comm_get_parent(MPI_Comm *parent);
int MPI_Comm_join(int fd, MPI_Comm *intercomm);
int MPI_Comm_spawn(const char *command, char *argv[], int maxprocs, MPI_Info info, int root,
                   MPI_Comm comm, MPI_Comm *intercomm, int array_of_errcodes[]);
int MPI_Comm_spawn_multiple(int count, char *array_of_commands[], char **array_of_argv[],
                            const int array_of_maxprocs[], const MPI_Info array_of_info[],
                            int root, MPI_Comm comm, MPI_Comm *intercomm, int array_of_errcodes[]);
int MPI_Lookup_name(const char *service_name, MPI_Info info, char *port_name);
int MPI_Open_port(MPI_Info info, char *port_name);
int MPI_Publish_name(const char *service_name, MPI_Info info, const char *port_name);
int MPI_Unpublish_name(const char *service_name, MPI_Info info, const char *port_name);
int MPI_Comm_set_info(MPI_Comm comm, MPI_Info info);
int MPI_Comm_get_info(MPI_Comm comm, MPI_Info *info);


int MPI_Accumulate(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
                   int target_rank, MPI_Aint target_disp, int target_count,
                   MPI_Datatype target_datatype, MPI_Op op, MPI_Win win)
                   ;
int MPI_Get(void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
            int target_rank, MPI_Aint target_disp, int target_count,
            MPI_Datatype target_datatype, MPI_Win win) ;
int MPI_Put(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
            int target_rank, MPI_Aint target_disp, int target_count,
            MPI_Datatype target_datatype, MPI_Win win) ;
int MPI_Win_complete(MPI_Win win);
int MPI_Win_create(void *base, MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm,
                   MPI_Win *win);
int MPI_Win_fence(int assert, MPI_Win win);
int MPI_Win_free(MPI_Win *win);
int MPI_Win_get_group(MPI_Win win, MPI_Group *group);
int MPI_Win_lock(int lock_type, int rank, int assert, MPI_Win win);
int MPI_Win_post(MPI_Group group, int assert, MPI_Win win);
int MPI_Win_start(MPI_Group group, int assert, MPI_Win win);
int MPI_Win_test(MPI_Win win, int *flag);
int MPI_Win_unlock(int rank, MPI_Win win);
int MPI_Win_wait(MPI_Win win);


int MPI_Win_allocate(MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm, void *baseptr,
                     MPI_Win *win);
int MPI_Win_allocate_shared(MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm,
                            void *baseptr, MPI_Win *win);
int MPI_Win_shared_query(MPI_Win win, int rank, MPI_Aint *size, int *disp_unit, void *baseptr);
int MPI_Win_create_dynamic(MPI_Info info, MPI_Comm comm, MPI_Win *win);
int MPI_Win_attach(MPI_Win win, void *base, MPI_Aint size);
int MPI_Win_detach(MPI_Win win, const void *base);
int MPI_Win_get_info(MPI_Win win, MPI_Info *info_used);
int MPI_Win_set_info(MPI_Win win, MPI_Info info);
int MPI_Get_accumulate(const void *origin_addr, int origin_count,
                        MPI_Datatype origin_datatype, void *result_addr, int result_count,
                        MPI_Datatype result_datatype, int target_rank, MPI_Aint target_disp,
                        int target_count, MPI_Datatype target_datatype, MPI_Op op, MPI_Win win)
                       
                        ;
int MPI_Fetch_and_op(const void *origin_addr, void *result_addr,
                      MPI_Datatype datatype, int target_rank, MPI_Aint target_disp,
                      MPI_Op op, MPI_Win win)
                      ;
int MPI_Compare_and_swap(const void *origin_addr, const void *compare_addr,
                          void *result_addr, MPI_Datatype datatype, int target_rank,
                          MPI_Aint target_disp, MPI_Win win)
                         
                         
                          ;
int MPI_Rput(const void *origin_addr, int origin_count,
              MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
              int target_count, MPI_Datatype target_datatype, MPI_Win win,
              MPI_Request *request)
              ;
int MPI_Rget(void *origin_addr, int origin_count,
              MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
              int target_count, MPI_Datatype target_datatype, MPI_Win win,
              MPI_Request *request)
              ;
int MPI_Raccumulate(const void *origin_addr, int origin_count,
                     MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
                     int target_count, MPI_Datatype target_datatype, MPI_Op op, MPI_Win win,
                     MPI_Request *request)
                     ;
int MPI_Rget_accumulate(const void *origin_addr, int origin_count,
                         MPI_Datatype origin_datatype, void *result_addr, int result_count,
                         MPI_Datatype result_datatype, int target_rank, MPI_Aint target_disp,
                         int target_count, MPI_Datatype target_datatype, MPI_Op op, MPI_Win win,
                         MPI_Request *request)
                        
                         ;
int MPI_Win_lock_all(int assert, MPI_Win win);
int MPI_Win_unlock_all(MPI_Win win);
int MPI_Win_flush(int rank, MPI_Win win);
int MPI_Win_flush_all(MPI_Win win);
int MPI_Win_flush_local(int rank, MPI_Win win);
int MPI_Win_flush_local_all(MPI_Win win);
int MPI_Win_sync(MPI_Win win);


int MPI_Add_error_class(int *errorclass);
int MPI_Add_error_code(int errorclass, int *errorcode);
int MPI_Add_error_string(int errorcode, const char *string);
int MPI_Comm_call_errhandler(MPI_Comm comm, int errorcode);
int MPI_Comm_create_keyval(MPI_Comm_copy_attr_function *comm_copy_attr_fn,
                           MPI_Comm_delete_attr_function *comm_delete_attr_fn, int *comm_keyval,
                           void *extra_state);
int MPI_Comm_delete_attr(MPI_Comm comm, int comm_keyval);
int MPI_Comm_free_keyval(int *comm_keyval);
int MPI_Comm_get_attr(MPI_Comm comm, int comm_keyval, void *attribute_val, int *flag);
int MPI_Comm_get_name(MPI_Comm comm, char *comm_name, int *resultlen);
int MPI_Comm_set_attr(MPI_Comm comm, int comm_keyval, void *attribute_val);
int MPI_Comm_set_name(MPI_Comm comm, const char *comm_name);
int MPI_File_call_errhandler(MPI_File fh, int errorcode);
int MPI_Grequest_complete(MPI_Request request);
int MPI_Grequest_start(MPI_Grequest_query_function *query_fn, MPI_Grequest_free_function *free_fn,
                       MPI_Grequest_cancel_function *cancel_fn, void *extra_state,
                       MPI_Request *request);
int MPI_Init_thread(int *argc, char ***argv, int required, int *provided);
int MPI_Is_thread_main(int *flag);
int MPI_Query_thread(int *provided);
int MPI_Status_set_cancelled(MPI_Status *status, int flag);
int MPI_Status_set_elements(MPI_Status *status, MPI_Datatype datatype, int count);
int MPI_Type_create_keyval(MPI_Type_copy_attr_function *type_copy_attr_fn,
                           MPI_Type_delete_attr_function *type_delete_attr_fn,
                           int *type_keyval, void *extra_state);
int MPI_Type_delete_attr(MPI_Datatype datatype, int type_keyval);
int MPI_Type_dup(MPI_Datatype oldtype, MPI_Datatype *newtype);
int MPI_Type_free_keyval(int *type_keyval);
int MPI_Type_get_attr(MPI_Datatype datatype, int type_keyval, void *attribute_val, int *flag);
int MPI_Type_get_contents(MPI_Datatype datatype, int max_integers, int max_addresses,
                          int max_datatypes, int array_of_integers[],
                          MPI_Aint array_of_addresses[], MPI_Datatype array_of_datatypes[]);
int MPI_Type_get_envelope(MPI_Datatype datatype, int *num_integers, int *num_addresses,
                          int *num_datatypes, int *combiner);
int MPI_Type_get_name(MPI_Datatype datatype, char *type_name, int *resultlen);
int MPI_Type_set_attr(MPI_Datatype datatype, int type_keyval, void *attribute_val);
int MPI_Type_set_name(MPI_Datatype datatype, const char *type_name);
int MPI_Type_match_size(int typeclass, int size, MPI_Datatype *datatype);
int MPI_Win_call_errhandler(MPI_Win win, int errorcode);
int MPI_Win_create_keyval(MPI_Win_copy_attr_function *win_copy_attr_fn,
                          MPI_Win_delete_attr_function *win_delete_attr_fn, int *win_keyval,
                          void *extra_state);
int MPI_Win_delete_attr(MPI_Win win, int win_keyval);
int MPI_Win_free_keyval(int *win_keyval);
int MPI_Win_get_attr(MPI_Win win, int win_keyval, void *attribute_val, int *flag);
int MPI_Win_get_name(MPI_Win win, char *win_name, int *resultlen);
int MPI_Win_set_attr(MPI_Win win, int win_keyval, void *attribute_val);
int MPI_Win_set_name(MPI_Win win, const char *win_name);

int MPI_Alloc_mem(MPI_Aint size, MPI_Info info, void *baseptr);
int MPI_Comm_create_errhandler(MPI_Comm_errhandler_function *comm_errhandler_fn,
                               MPI_Errhandler *errhandler);
int MPI_Comm_get_errhandler(MPI_Comm comm, MPI_Errhandler *errhandler);
int MPI_Comm_set_errhandler(MPI_Comm comm, MPI_Errhandler errhandler);
int MPI_File_create_errhandler(MPI_File_errhandler_function *file_errhandler_fn,
                               MPI_Errhandler *errhandler);
int MPI_File_get_errhandler(MPI_File file, MPI_Errhandler *errhandler);
int MPI_File_set_errhandler(MPI_File file, MPI_Errhandler errhandler);
int MPI_Finalized(int *flag);
int MPI_Free_mem(void *base);
int MPI_Get_address(const void *location, MPI_Aint *address);
int MPI_Info_create(MPI_Info *info);
int MPI_Info_delete(MPI_Info info, const char *key);
int MPI_Info_dup(MPI_Info info, MPI_Info *newinfo);
int MPI_Info_free(MPI_Info *info);
int MPI_Info_get(MPI_Info info, const char *key, int valuelen, char *value, int *flag);
int MPI_Info_get_nkeys(MPI_Info info, int *nkeys);
int MPI_Info_get_nthkey(MPI_Info info, int n, char *key);
int MPI_Info_get_valuelen(MPI_Info info, const char *key, int *valuelen, int *flag);
int MPI_Info_set(MPI_Info info, const char *key, const char *value);
int MPI_Pack_external(const char datarep[], const void *inbuf, int incount,
                      MPI_Datatype datatype, void *outbuf, MPI_Aint outsize, MPI_Aint *position)
                      ;
int MPI_Pack_external_size(const char datarep[], int incount, MPI_Datatype datatype,
                           MPI_Aint *size);
int MPI_Request_get_status(MPI_Request request, int *flag, MPI_Status *status);
int MPI_Status_c2f(const MPI_Status *c_status, MPI_Fint *f_status);
int MPI_Status_f2c(const MPI_Fint *f_status, MPI_Status *c_status);
int MPI_Type_create_darray(int size, int rank, int ndims, const int array_of_gsizes[],
                           const int array_of_distribs[], const int array_of_dargs[],
                           const int array_of_psizes[], int order, MPI_Datatype oldtype,
                           MPI_Datatype *newtype);
int MPI_Type_create_hindexed(int count, const int array_of_blocklengths[],
                             const MPI_Aint array_of_displacements[], MPI_Datatype oldtype,
                             MPI_Datatype *newtype);
int MPI_Type_create_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype,
                            MPI_Datatype *newtype);
int MPI_Type_create_indexed_block(int count, int blocklength, const int array_of_displacements[],
                                  MPI_Datatype oldtype, MPI_Datatype *newtype);
int MPI_Type_create_hindexed_block(int count, int blocklength,
                                   const MPI_Aint array_of_displacements[],
                                   MPI_Datatype oldtype, MPI_Datatype *newtype);
int MPI_Type_create_resized(MPI_Datatype oldtype, MPI_Aint lb, MPI_Aint extent,
                            MPI_Datatype *newtype);
int MPI_Type_create_struct(int count, const int array_of_blocklengths[],
                           const MPI_Aint array_of_displacements[],
                           const MPI_Datatype array_of_types[], MPI_Datatype *newtype);
int MPI_Type_create_subarray(int ndims, const int array_of_sizes[],
                             const int array_of_subsizes[], const int array_of_starts[],
                             int order, MPI_Datatype oldtype, MPI_Datatype *newtype);
int MPI_Type_get_extent(MPI_Datatype datatype, MPI_Aint *lb, MPI_Aint *extent);
int MPI_Type_get_true_extent(MPI_Datatype datatype, MPI_Aint *true_lb, MPI_Aint *true_extent);
int MPI_Unpack_external(const char datarep[], const void *inbuf, MPI_Aint insize,
                        MPI_Aint *position, void *outbuf, int outcount, MPI_Datatype datatype)
                        ;
int MPI_Win_create_errhandler(MPI_Win_errhandler_function *win_errhandler_fn,
                              MPI_Errhandler *errhandler);
int MPI_Win_get_errhandler(MPI_Win win, MPI_Errhandler *errhandler);
int MPI_Win_set_errhandler(MPI_Win win, MPI_Errhandler errhandler);




int MPI_Type_create_f90_integer(int range, MPI_Datatype *newtype);
int MPI_Type_create_f90_real(int precision, int range, MPI_Datatype *newtype);
int MPI_Type_create_f90_complex(int precision, int range, MPI_Datatype *newtype);

int MPI_Reduce_local(const void *inbuf, void *inoutbuf, int count, MPI_Datatype datatype,
                     MPI_Op op)
                     ;
int MPI_Op_commutative(MPI_Op op, int *commute);
int MPI_Reduce_scatter_block(const void *sendbuf, void *recvbuf, int recvcount,
                             MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                            
                             ;
int MPI_Dist_graph_create_adjacent(MPI_Comm comm_old, int indegree, const int sources[],
                                   const int sourceweights[], int outdegree,
                                   const int destinations[], const int destweights[],
                                   MPI_Info info, int reorder, MPI_Comm *comm_dist_graph);
int MPI_Dist_graph_create(MPI_Comm comm_old, int n, const int sources[], const int degrees[],
                          const int destinations[], const int weights[], MPI_Info info,
                          int reorder, MPI_Comm *comm_dist_graph);
int MPI_Dist_graph_neighbors_count(MPI_Comm comm, int *indegree, int *outdegree, int *weighted);
int MPI_Dist_graph_neighbors(MPI_Comm comm, int maxindegree, int sources[], int sourceweights[],
                             int maxoutdegree, int destinations[], int destweights[]);


int MPI_Improbe(int source, int tag, MPI_Comm comm, int *flag, MPI_Message *message,
                MPI_Status *status);
int MPI_Imrecv(void *buf, int count, MPI_Datatype datatype, MPI_Message *message,
               MPI_Request *request) ;
int MPI_Mprobe(int source, int tag, MPI_Comm comm, MPI_Message *message, MPI_Status *status);
int MPI_Mrecv(void *buf, int count, MPI_Datatype datatype, MPI_Message *message,
              MPI_Status *status) ;


int MPI_Comm_idup(MPI_Comm comm, MPI_Comm *newcomm, MPI_Request *request);
int MPI_Ibarrier(MPI_Comm comm, MPI_Request *request);
int MPI_Ibcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm,
               MPI_Request *request) ;
int MPI_Igather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm,
                MPI_Request *request)
                ;
int MPI_Igatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                 const int recvcounts[], const int displs[], MPI_Datatype recvtype, int root,
                 MPI_Comm comm, MPI_Request *request)
                 ;
int MPI_Iscatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                 int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm,
                 MPI_Request *request)
                 ;
int MPI_Iscatterv(const void *sendbuf, const int sendcounts[], const int displs[],
                  MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                  int root, MPI_Comm comm, MPI_Request *request)
                  ;
int MPI_Iallgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                   int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                   ;
int MPI_Iallgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                    const int recvcounts[], const int displs[], MPI_Datatype recvtype,
                    MPI_Comm comm, MPI_Request *request)
                    ;
int MPI_Ialltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                  int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                  ;
int MPI_Ialltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[],
                   MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                   const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm,
                   MPI_Request *request)
                   ;
int MPI_Ialltoallw(const void *sendbuf, const int sendcounts[], const int sdispls[],
                   const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
                   const int rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm,
                   MPI_Request *request);
int MPI_Ireduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                MPI_Op op, int root, MPI_Comm comm, MPI_Request *request)
                ;
int MPI_Iallreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                   MPI_Op op, MPI_Comm comm, MPI_Request *request)
                   ;
int MPI_Ireduce_scatter(const void *sendbuf, void *recvbuf, const int recvcounts[],
                        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request)
                        ;
int MPI_Ireduce_scatter_block(const void *sendbuf, void *recvbuf, int recvcount,
                              MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                              MPI_Request *request)
                             
                              ;
int MPI_Iscan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
              MPI_Comm comm, MPI_Request *request)
              ;
int MPI_Iexscan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                MPI_Op op, MPI_Comm comm, MPI_Request *request)
                ;


int MPI_Ineighbor_allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                            void *recvbuf, int recvcount, MPI_Datatype recvtype,
                            MPI_Comm comm, MPI_Request *request)
                           
                            ;
int MPI_Ineighbor_allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                             void *recvbuf, const int recvcounts[], const int displs[],
                             MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                            
                             ;
int MPI_Ineighbor_alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                           void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                           MPI_Request *request)
                          
                           ;
int MPI_Ineighbor_alltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[],
                            MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                            const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm,
                            MPI_Request *request)
                           
                            ;
int MPI_Ineighbor_alltoallw(const void *sendbuf, const int sendcounts[],
                            const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
                            void *recvbuf, const int recvcounts[], const MPI_Aint rdispls[],
                            const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Request *request);
int MPI_Neighbor_allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                           void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                          
                           ;
int MPI_Neighbor_allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                            void *recvbuf, const int recvcounts[], const int displs[],
                            MPI_Datatype recvtype, MPI_Comm comm)
                           
                            ;
int MPI_Neighbor_alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                          void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                         
                          ;
int MPI_Neighbor_alltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[],
                           MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                           const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm)
                          
                           ;
int MPI_Neighbor_alltoallw(const void *sendbuf, const int sendcounts[], const MPI_Aint sdispls[],
                           const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
                           const MPI_Aint rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm);


int MPI_Comm_split_type(MPI_Comm comm, int split_type, int key, MPI_Info info, MPI_Comm *newcomm);


int MPI_Get_elements_x(const MPI_Status *status, MPI_Datatype datatype, MPI_Count *count);
int MPI_Status_set_elements_x(MPI_Status *status, MPI_Datatype datatype, MPI_Count count);
int MPI_Type_get_extent_x(MPI_Datatype datatype, MPI_Count *lb, MPI_Count *extent);
int MPI_Type_get_true_extent_x(MPI_Datatype datatype, MPI_Count *lb, MPI_Count *extent);
int MPI_Type_size_x(MPI_Datatype datatype, MPI_Count *size);


int MPI_Comm_create_group(MPI_Comm comm, MPI_Group group, int tag, MPI_Comm *newcomm);


MPI_Aint MPI_Aint_add(MPI_Aint base, MPI_Aint disp);
MPI_Aint MPI_Aint_diff(MPI_Aint addr1, MPI_Aint addr2);





int MPI_T_init_thread(int required, int *provided);
int MPI_T_finalize(void);
int MPI_T_enum_get_info(MPI_T_enum enumtype, int *num, char *name, int *name_len);
int MPI_T_enum_get_item(MPI_T_enum enumtype, int indx, int *value, char *name, int *name_len);
int MPI_T_cvar_get_num(int *num_cvar);
int MPI_T_cvar_get_info(int cvar_index, char *name, int *name_len, int *verbosity,
                        MPI_Datatype *datatype, MPI_T_enum *enumtype, char *desc, int *desc_len,
                        int *binding, int *scope);
int MPI_T_cvar_handle_alloc(int cvar_index, void *obj_handle, MPI_T_cvar_handle *handle,
                            int *count);
int MPI_T_cvar_handle_free(MPI_T_cvar_handle *handle);
int MPI_T_cvar_read(MPI_T_cvar_handle handle, void *buf);
int MPI_T_cvar_write(MPI_T_cvar_handle handle, const void *buf);
int MPI_T_pvar_get_num(int *num_pvar);
int MPI_T_pvar_get_info(int pvar_index, char *name, int *name_len, int *verbosity, int *var_class,
                        MPI_Datatype *datatype, MPI_T_enum *enumtype, char *desc, int *desc_len,
                        int *binding, int *readonly, int *continuous, int *atomic);
int MPI_T_pvar_session_create(MPI_T_pvar_session *session);
int MPI_T_pvar_session_free(MPI_T_pvar_session *session);
int MPI_T_pvar_handle_alloc(MPI_T_pvar_session session, int pvar_index, void *obj_handle,
                            MPI_T_pvar_handle *handle, int *count);
int MPI_T_pvar_handle_free(MPI_T_pvar_session session, MPI_T_pvar_handle *handle);
int MPI_T_pvar_start(MPI_T_pvar_session session, MPI_T_pvar_handle handle);
int MPI_T_pvar_stop(MPI_T_pvar_session session, MPI_T_pvar_handle handle);
int MPI_T_pvar_read(MPI_T_pvar_session session, MPI_T_pvar_handle handle, void *buf);
int MPI_T_pvar_write(MPI_T_pvar_session session, MPI_T_pvar_handle handle, const void *buf);
int MPI_T_pvar_reset(MPI_T_pvar_session session, MPI_T_pvar_handle handle);
int MPI_T_pvar_readreset(MPI_T_pvar_session session, MPI_T_pvar_handle handle, void *buf);
int MPI_T_category_get_num(int *num_cat);
int MPI_T_category_get_info(int cat_index, char *name, int *name_len, char *desc, int *desc_len,
                            int *num_cvars, int *num_pvars, int *num_categories);
int MPI_T_category_get_cvars(int cat_index, int len, int indices[]);
int MPI_T_category_get_pvars(int cat_index, int len, int indices[]);
int MPI_T_category_get_categories(int cat_index, int len, int indices[]);
int MPI_T_category_changed(int *stamp);
int MPI_T_cvar_get_index(const char *name, int *cvar_index);
int MPI_T_pvar_get_index(const char *name, int var_class, int *pvar_index);
int MPI_T_category_get_index(const char *name, int *cat_index);





int MPIX_Comm_failure_ack(MPI_Comm comm);
int MPIX_Comm_failure_get_acked(MPI_Comm comm, MPI_Group *failedgrp);
int MPIX_Comm_revoke(MPI_Comm comm);
int MPIX_Comm_shrink(MPI_Comm comm, MPI_Comm *newcomm);
int MPIX_Comm_agree(MPI_Comm comm, int *flag);
# 1661 "/sw/eb/sw/impi/2018.5.288-iccifort-2019.5.281/include64/mpi.h"
int PMPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm) ;
int PMPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
              MPI_Comm comm, MPI_Status *status) ;
int PMPI_Get_count(const MPI_Status *status, MPI_Datatype datatype, int *count);
int PMPI_Bsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
               MPI_Comm comm) ;
int PMPI_Ssend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
               MPI_Comm comm) ;
int PMPI_Rsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
               MPI_Comm comm) ;
int PMPI_Buffer_attach(void *buffer, int size);
int PMPI_Buffer_detach(void *buffer_addr, int *size);
int PMPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
               MPI_Comm comm, MPI_Request *request) ;
int PMPI_Ibsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                MPI_Comm comm, MPI_Request *request) ;
int PMPI_Issend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                MPI_Comm comm, MPI_Request *request) ;
int PMPI_Irsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                MPI_Comm comm, MPI_Request *request) ;
int PMPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
               MPI_Comm comm, MPI_Request *request) ;
int PMPI_Wait(MPI_Request *request, MPI_Status *status);
int PMPI_Test(MPI_Request *request, int *flag, MPI_Status *status);
int PMPI_Request_free(MPI_Request *request);
int PMPI_Waitany(int count, MPI_Request array_of_requests[], int *indx, MPI_Status *status);
int PMPI_Testany(int count, MPI_Request array_of_requests[], int *indx, int *flag,
                 MPI_Status *status);
int PMPI_Waitall(int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]);
int PMPI_Testall(int count, MPI_Request array_of_requests[], int *flag,
                 MPI_Status array_of_statuses[]);
int PMPI_Waitsome(int incount, MPI_Request array_of_requests[], int *outcount,
                  int array_of_indices[], MPI_Status array_of_statuses[]);
int PMPI_Testsome(int incount, MPI_Request array_of_requests[], int *outcount,
                  int array_of_indices[], MPI_Status array_of_statuses[]);
int PMPI_Iprobe(int source, int tag, MPI_Comm comm, int *flag, MPI_Status *status);
int PMPI_Probe(int source, int tag, MPI_Comm comm, MPI_Status *status);
int PMPI_Cancel(MPI_Request *request);
int PMPI_Test_cancelled(const MPI_Status *status, int *flag);
int PMPI_Send_init(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                   MPI_Comm comm, MPI_Request *request) ;
int PMPI_Bsend_init(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                    MPI_Comm comm, MPI_Request *request) ;
int PMPI_Ssend_init(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                    MPI_Comm comm, MPI_Request *request) ;
int PMPI_Rsend_init(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                    MPI_Comm comm, MPI_Request *request) ;
int PMPI_Recv_init(void *buf, int count, MPI_Datatype datatype, int source, int tag,
                   MPI_Comm comm, MPI_Request *request) ;
int PMPI_Start(MPI_Request *request);
int PMPI_Startall(int count, MPI_Request array_of_requests[]);
int PMPI_Sendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, int dest,
                  int sendtag, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                  int source, int recvtag, MPI_Comm comm, MPI_Status *status)
                 
                  ;
int PMPI_Sendrecv_replace(void *buf, int count, MPI_Datatype datatype, int dest,
                          int sendtag, int source, int recvtag, MPI_Comm comm,
                          MPI_Status *status) ;
int PMPI_Type_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype *newtype);
int PMPI_Type_vector(int count, int blocklength, int stride, MPI_Datatype oldtype,
                     MPI_Datatype *newtype);
int PMPI_Type_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype,
                      MPI_Datatype *newtype);
int PMPI_Type_indexed(int count, const int *array_of_blocklengths,
                      const int *array_of_displacements, MPI_Datatype oldtype,
                      MPI_Datatype *newtype);
int PMPI_Type_hindexed(int count, const int *array_of_blocklengths,
                       const MPI_Aint *array_of_displacements, MPI_Datatype oldtype,
                       MPI_Datatype *newtype);
int PMPI_Type_struct(int count, const int *array_of_blocklengths,
                     const MPI_Aint *array_of_displacements,
                     const MPI_Datatype *array_of_types, MPI_Datatype *newtype);
int PMPI_Address(const void *location, MPI_Aint *address);
int PMPI_Type_extent(MPI_Datatype datatype, MPI_Aint *extent);
int PMPI_Type_size(MPI_Datatype datatype, int *size);
int PMPI_Type_lb(MPI_Datatype datatype, MPI_Aint *displacement);
int PMPI_Type_ub(MPI_Datatype datatype, MPI_Aint *displacement);
int PMPI_Type_commit(MPI_Datatype *datatype);
int PMPI_Type_free(MPI_Datatype *datatype);
int PMPI_Get_elements(const MPI_Status *status, MPI_Datatype datatype, int *count);
int PMPI_Pack(const void *inbuf, int incount, MPI_Datatype datatype, void *outbuf,
              int outsize, int *position, MPI_Comm comm) ;
int PMPI_Unpack(const void *inbuf, int insize, int *position, void *outbuf, int outcount,
                MPI_Datatype datatype, MPI_Comm comm) ;
int PMPI_Pack_size(int incount, MPI_Datatype datatype, MPI_Comm comm, int *size);
int PMPI_Barrier(MPI_Comm comm);
int PMPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)
               ;
int PMPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
                ;
int PMPI_Gatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                 const int *recvcounts, const int *displs, MPI_Datatype recvtype, int root,
                 MPI_Comm comm)
                 ;
int PMPI_Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                 int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
                 ;
int PMPI_Scatterv(const void *sendbuf, const int *sendcounts, const int *displs,
                  MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                  int root, MPI_Comm comm)
                  ;
int PMPI_Allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                   int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                   ;
int PMPI_Allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                    const int *recvcounts, const int *displs, MPI_Datatype recvtype, MPI_Comm comm)
                    ;
int PMPI_Alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                  int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                  ;
int PMPI_Alltoallv(const void *sendbuf, const int *sendcounts, const int *sdispls,
                   MPI_Datatype sendtype, void *recvbuf, const int *recvcounts,
                   const int *rdispls, MPI_Datatype recvtype, MPI_Comm comm)
                   ;
int PMPI_Alltoallw(const void *sendbuf, const int sendcounts[], const int sdispls[],
                   const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
                   const int rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm);
int PMPI_Exscan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                MPI_Op op, MPI_Comm comm)
                ;
int PMPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                MPI_Op op, int root, MPI_Comm comm)
                ;
int PMPI_Op_create(MPI_User_function *user_fn, int commute, MPI_Op *op);
int PMPI_Op_free(MPI_Op *op);
int PMPI_Allreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                   MPI_Op op, MPI_Comm comm)
                   ;
int PMPI_Reduce_scatter(const void *sendbuf, void *recvbuf, const int recvcounts[],
                        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                        ;
int PMPI_Scan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
              MPI_Comm comm)
              ;
int PMPI_Group_size(MPI_Group group, int *size);
int PMPI_Group_rank(MPI_Group group, int *rank);
int PMPI_Group_translate_ranks(MPI_Group group1, int n, const int ranks1[], MPI_Group group2,
                               int ranks2[]);
int PMPI_Group_compare(MPI_Group group1, MPI_Group group2, int *result);
int PMPI_Comm_group(MPI_Comm comm, MPI_Group *group);
int PMPI_Group_union(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup);
int PMPI_Group_intersection(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup);
int PMPI_Group_difference(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup);
int PMPI_Group_incl(MPI_Group group, int n, const int ranks[], MPI_Group *newgroup);
int PMPI_Group_excl(MPI_Group group, int n, const int ranks[], MPI_Group *newgroup);
int PMPI_Group_range_incl(MPI_Group group, int n, int ranges[][3], MPI_Group *newgroup);
int PMPI_Group_range_excl(MPI_Group group, int n, int ranges[][3], MPI_Group *newgroup);
int PMPI_Group_free(MPI_Group *group);
int PMPI_Comm_size(MPI_Comm comm, int *size);
int PMPI_Comm_rank(MPI_Comm comm, int *rank);
int PMPI_Comm_compare(MPI_Comm comm1, MPI_Comm comm2, int *result);
int PMPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm);
int PMPI_Comm_dup_with_info(MPI_Comm comm, MPI_Info info, MPI_Comm *newcomm);
int PMPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm);
int PMPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm);
int PMPI_Comm_free(MPI_Comm *comm);
int PMPI_Comm_test_inter(MPI_Comm comm, int *flag);
int PMPI_Comm_remote_size(MPI_Comm comm, int *size);
int PMPI_Comm_remote_group(MPI_Comm comm, MPI_Group *group);
int PMPI_Intercomm_create(MPI_Comm local_comm, int local_leader, MPI_Comm peer_comm,
                          int remote_leader, int tag, MPI_Comm *newintercomm);
int PMPI_Intercomm_merge(MPI_Comm intercomm, int high, MPI_Comm *newintracomm);
int PMPI_Keyval_create(MPI_Copy_function *copy_fn, MPI_Delete_function *delete_fn,
                       int *keyval, void *extra_state);
int PMPI_Keyval_free(int *keyval);
int PMPI_Attr_put(MPI_Comm comm, int keyval, void *attribute_val);
int PMPI_Attr_get(MPI_Comm comm, int keyval, void *attribute_val, int *flag);
int PMPI_Attr_delete(MPI_Comm comm, int keyval);
int PMPI_Topo_test(MPI_Comm comm, int *status);
int PMPI_Cart_create(MPI_Comm comm_old, int ndims, const int dims[], const int periods[],
                     int reorder, MPI_Comm *comm_cart);
int PMPI_Dims_create(int nnodes, int ndims, int dims[]);
int PMPI_Graph_create(MPI_Comm comm_old, int nnodes, const int indx[], const int edges[],
                      int reorder, MPI_Comm *comm_graph);
int PMPI_Graphdims_get(MPI_Comm comm, int *nnodes, int *nedges);
int PMPI_Graph_get(MPI_Comm comm, int maxindex, int maxedges, int indx[], int edges[]);
int PMPI_Cartdim_get(MPI_Comm comm, int *ndims);
int PMPI_Cart_get(MPI_Comm comm, int maxdims, int dims[], int periods[], int coords[]);
int PMPI_Cart_rank(MPI_Comm comm, const int coords[], int *rank);
int PMPI_Cart_coords(MPI_Comm comm, int rank, int maxdims, int coords[]);
int PMPI_Graph_neighbors_count(MPI_Comm comm, int rank, int *nneighbors);
int PMPI_Graph_neighbors(MPI_Comm comm, int rank, int maxneighbors, int neighbors[]);
int PMPI_Cart_shift(MPI_Comm comm, int direction, int disp, int *rank_source, int *rank_dest);
int PMPI_Cart_sub(MPI_Comm comm, const int remain_dims[], MPI_Comm *newcomm);
int PMPI_Cart_map(MPI_Comm comm, int ndims, const int dims[], const int periods[], int *newrank);
int PMPI_Graph_map(MPI_Comm comm, int nnodes, const int indx[], const int edges[], int *newrank);
int PMPI_Get_processor_name(char *name, int *resultlen);
int PMPI_Get_version(int *version, int *subversion);
int PMPI_Get_library_version(char *version, int *resultlen);
int PMPI_Errhandler_create(MPI_Handler_function *function, MPI_Errhandler *errhandler);
int PMPI_Errhandler_set(MPI_Comm comm, MPI_Errhandler errhandler);
int PMPI_Errhandler_get(MPI_Comm comm, MPI_Errhandler *errhandler);
int PMPI_Errhandler_free(MPI_Errhandler *errhandler);
int PMPI_Error_string(int errorcode, char *string, int *resultlen);
int PMPI_Error_class(int errorcode, int *errorclass);
double PMPI_Wtime(void);
double PMPI_Wtick(void);
int PMPI_Init(int *argc, char ***argv);
int PMPI_Finalize(void);
int PMPI_Initialized(int *flag);
int PMPI_Abort(MPI_Comm comm, int errorcode);



int PMPI_Pcontrol(const int level, ...);


int PMPI_Close_port(const char *port_name);
int PMPI_Comm_accept(const char *port_name, MPI_Info info, int root, MPI_Comm comm,
                     MPI_Comm *newcomm);
int PMPI_Comm_connect(const char *port_name, MPI_Info info, int root, MPI_Comm comm,
                      MPI_Comm *newcomm);
int PMPI_Comm_disconnect(MPI_Comm *comm);
int PMPI_Comm_get_parent(MPI_Comm *parent);
int PMPI_Comm_join(int fd, MPI_Comm *intercomm);
int PMPI_Comm_spawn(const char *command, char *argv[], int maxprocs, MPI_Info info, int root,
                    MPI_Comm comm, MPI_Comm *intercomm, int array_of_errcodes[]);
int PMPI_Comm_spawn_multiple(int count, char *array_of_commands[], char **array_of_argv[],
                             const int array_of_maxprocs[], const MPI_Info array_of_info[],
                             int root, MPI_Comm comm, MPI_Comm *intercomm, int array_of_errcodes[]);
int PMPI_Lookup_name(const char *service_name, MPI_Info info, char *port_name);
int PMPI_Open_port(MPI_Info info, char *port_name);
int PMPI_Publish_name(const char *service_name, MPI_Info info, const char *port_name);
int PMPI_Unpublish_name(const char *service_name, MPI_Info info, const char *port_name);
int PMPI_Comm_set_info(MPI_Comm comm, MPI_Info info);
int PMPI_Comm_get_info(MPI_Comm comm, MPI_Info *info);


int PMPI_Accumulate(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
                    int target_rank, MPI_Aint target_disp, int target_count,
                    MPI_Datatype target_datatype, MPI_Op op, MPI_Win win)
                    ;
int PMPI_Get(void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
             int target_rank, MPI_Aint target_disp, int target_count,
             MPI_Datatype target_datatype, MPI_Win win) ;
int PMPI_Put(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
             int target_rank, MPI_Aint target_disp, int target_count,
             MPI_Datatype target_datatype, MPI_Win win) ;
int PMPI_Win_complete(MPI_Win win);
int PMPI_Win_create(void *base, MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm,
                    MPI_Win *win);
int PMPI_Win_fence(int assert, MPI_Win win);
int PMPI_Win_free(MPI_Win *win);
int PMPI_Win_get_group(MPI_Win win, MPI_Group *group);
int PMPI_Win_lock(int lock_type, int rank, int assert, MPI_Win win);
int PMPI_Win_post(MPI_Group group, int assert, MPI_Win win);
int PMPI_Win_start(MPI_Group group, int assert, MPI_Win win);
int PMPI_Win_test(MPI_Win win, int *flag);
int PMPI_Win_unlock(int rank, MPI_Win win);
int PMPI_Win_wait(MPI_Win win);


int PMPI_Win_allocate(MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm, void *baseptr,
                      MPI_Win *win);
int PMPI_Win_allocate_shared(MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm,
                             void *baseptr, MPI_Win *win);
int PMPI_Win_shared_query(MPI_Win win, int rank, MPI_Aint *size, int *disp_unit, void *baseptr);
int PMPI_Win_create_dynamic(MPI_Info info, MPI_Comm comm, MPI_Win *win);
int PMPI_Win_attach(MPI_Win win, void *base, MPI_Aint size);
int PMPI_Win_detach(MPI_Win win, const void *base);
int PMPI_Win_get_info(MPI_Win win, MPI_Info *info_used);
int PMPI_Win_set_info(MPI_Win win, MPI_Info info);
int PMPI_Get_accumulate(const void *origin_addr, int origin_count,
                         MPI_Datatype origin_datatype, void *result_addr, int result_count,
                         MPI_Datatype result_datatype, int target_rank, MPI_Aint target_disp,
                         int target_count, MPI_Datatype target_datatype, MPI_Op op, MPI_Win win)
                        
                         ;
int PMPI_Fetch_and_op(const void *origin_addr, void *result_addr,
                       MPI_Datatype datatype, int target_rank, MPI_Aint target_disp,
                       MPI_Op op, MPI_Win win)
                       ;
int PMPI_Compare_and_swap(const void *origin_addr, const void *compare_addr,
                           void *result_addr, MPI_Datatype datatype, int target_rank,
                           MPI_Aint target_disp, MPI_Win win)
                          
                          
                           ;
int PMPI_Rput(const void *origin_addr, int origin_count,
               MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
               int target_count, MPI_Datatype target_datatype, MPI_Win win,
               MPI_Request *request)
               ;
int PMPI_Rget(void *origin_addr, int origin_count,
               MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
               int target_count, MPI_Datatype target_datatype, MPI_Win win,
               MPI_Request *request)
               ;
int PMPI_Raccumulate(const void *origin_addr, int origin_count,
                      MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
                      int target_count, MPI_Datatype target_datatype, MPI_Op op, MPI_Win win,
                      MPI_Request *request)
                      ;
int PMPI_Rget_accumulate(const void *origin_addr, int origin_count,
                          MPI_Datatype origin_datatype, void *result_addr, int result_count,
                          MPI_Datatype result_datatype, int target_rank, MPI_Aint target_disp,
                          int target_count, MPI_Datatype target_datatype, MPI_Op op, MPI_Win win,
                          MPI_Request *request)
                         
                          ;
int PMPI_Win_lock_all(int assert, MPI_Win win);
int PMPI_Win_unlock_all(MPI_Win win);
int PMPI_Win_flush(int rank, MPI_Win win);
int PMPI_Win_flush_all(MPI_Win win);
int PMPI_Win_flush_local(int rank, MPI_Win win);
int PMPI_Win_flush_local_all(MPI_Win win);
int PMPI_Win_sync(MPI_Win win);


int PMPI_Add_error_class(int *errorclass);
int PMPI_Add_error_code(int errorclass, int *errorcode);
int PMPI_Add_error_string(int errorcode, const char *string);
int PMPI_Comm_call_errhandler(MPI_Comm comm, int errorcode);
int PMPI_Comm_create_keyval(MPI_Comm_copy_attr_function *comm_copy_attr_fn,
                            MPI_Comm_delete_attr_function *comm_delete_attr_fn, int *comm_keyval,
                            void *extra_state);
int PMPI_Comm_delete_attr(MPI_Comm comm, int comm_keyval);
int PMPI_Comm_free_keyval(int *comm_keyval);
int PMPI_Comm_get_attr(MPI_Comm comm, int comm_keyval, void *attribute_val, int *flag);
int PMPI_Comm_get_name(MPI_Comm comm, char *comm_name, int *resultlen);
int PMPI_Comm_set_attr(MPI_Comm comm, int comm_keyval, void *attribute_val);
int PMPI_Comm_set_name(MPI_Comm comm, const char *comm_name);
int PMPI_File_call_errhandler(MPI_File fh, int errorcode);
int PMPI_Grequest_complete(MPI_Request request);
int PMPI_Grequest_start(MPI_Grequest_query_function *query_fn, MPI_Grequest_free_function *free_fn,
                        MPI_Grequest_cancel_function *cancel_fn, void *extra_state,
                        MPI_Request *request);
int PMPI_Init_thread(int *argc, char ***argv, int required, int *provided);
int PMPI_Is_thread_main(int *flag);
int PMPI_Query_thread(int *provided);
int PMPI_Status_set_cancelled(MPI_Status *status, int flag);
int PMPI_Status_set_elements(MPI_Status *status, MPI_Datatype datatype, int count);
int PMPI_Type_create_keyval(MPI_Type_copy_attr_function *type_copy_attr_fn,
                            MPI_Type_delete_attr_function *type_delete_attr_fn,
                            int *type_keyval, void *extra_state);
int PMPI_Type_delete_attr(MPI_Datatype datatype, int type_keyval);
int PMPI_Type_dup(MPI_Datatype oldtype, MPI_Datatype *newtype);
int PMPI_Type_free_keyval(int *type_keyval);
int PMPI_Type_get_attr(MPI_Datatype datatype, int type_keyval, void *attribute_val, int *flag);
int PMPI_Type_get_contents(MPI_Datatype datatype, int max_integers, int max_addresses,
                           int max_datatypes, int array_of_integers[],
                           MPI_Aint array_of_addresses[], MPI_Datatype array_of_datatypes[]);
int PMPI_Type_get_envelope(MPI_Datatype datatype, int *num_integers, int *num_addresses,
                           int *num_datatypes, int *combiner);
int PMPI_Type_get_name(MPI_Datatype datatype, char *type_name, int *resultlen);
int PMPI_Type_set_attr(MPI_Datatype datatype, int type_keyval, void *attribute_val);
int PMPI_Type_set_name(MPI_Datatype datatype, const char *type_name);
int PMPI_Type_match_size(int typeclass, int size, MPI_Datatype *datatype);
int PMPI_Win_call_errhandler(MPI_Win win, int errorcode);
int PMPI_Win_create_keyval(MPI_Win_copy_attr_function *win_copy_attr_fn,
                           MPI_Win_delete_attr_function *win_delete_attr_fn, int *win_keyval,
                           void *extra_state);
int PMPI_Win_delete_attr(MPI_Win win, int win_keyval);
int PMPI_Win_free_keyval(int *win_keyval);
int PMPI_Win_get_attr(MPI_Win win, int win_keyval, void *attribute_val, int *flag);
int PMPI_Win_get_name(MPI_Win win, char *win_name, int *resultlen);
int PMPI_Win_set_attr(MPI_Win win, int win_keyval, void *attribute_val);
int PMPI_Win_set_name(MPI_Win win, const char *win_name);

int PMPI_Alloc_mem(MPI_Aint size, MPI_Info info, void *baseptr);
int PMPI_Comm_create_errhandler(MPI_Comm_errhandler_function *comm_errhandler_fn,
                                MPI_Errhandler *errhandler);
int PMPI_Comm_get_errhandler(MPI_Comm comm, MPI_Errhandler *errhandler);
int PMPI_Comm_set_errhandler(MPI_Comm comm, MPI_Errhandler errhandler);
int PMPI_File_create_errhandler(MPI_File_errhandler_function *file_errhandler_fn,
                                MPI_Errhandler *errhandler);
int PMPI_File_get_errhandler(MPI_File file, MPI_Errhandler *errhandler);
int PMPI_File_set_errhandler(MPI_File file, MPI_Errhandler errhandler);
int PMPI_Finalized(int *flag);
int PMPI_Free_mem(void *base);
int PMPI_Get_address(const void *location, MPI_Aint *address);
int PMPI_Info_create(MPI_Info *info);
int PMPI_Info_delete(MPI_Info info, const char *key);
int PMPI_Info_dup(MPI_Info info, MPI_Info *newinfo);
int PMPI_Info_free(MPI_Info *info);
int PMPI_Info_get(MPI_Info info, const char *key, int valuelen, char *value, int *flag);
int PMPI_Info_get_nkeys(MPI_Info info, int *nkeys);
int PMPI_Info_get_nthkey(MPI_Info info, int n, char *key);
int PMPI_Info_get_valuelen(MPI_Info info, const char *key, int *valuelen, int *flag);
int PMPI_Info_set(MPI_Info info, const char *key, const char *value);
int PMPI_Pack_external(const char datarep[], const void *inbuf, int incount,
                       MPI_Datatype datatype, void *outbuf, MPI_Aint outsize, MPI_Aint *position)
                       ;
int PMPI_Pack_external_size(const char datarep[], int incount, MPI_Datatype datatype,
                            MPI_Aint *size);
int PMPI_Request_get_status(MPI_Request request, int *flag, MPI_Status *status);
int PMPI_Status_c2f(const MPI_Status *c_status, MPI_Fint *f_status);
int PMPI_Status_f2c(const MPI_Fint *f_status, MPI_Status *c_status);
int PMPI_Type_create_darray(int size, int rank, int ndims, const int array_of_gsizes[],
                            const int array_of_distribs[], const int array_of_dargs[],
                            const int array_of_psizes[], int order, MPI_Datatype oldtype,
                            MPI_Datatype *newtype);
int PMPI_Type_create_hindexed(int count, const int array_of_blocklengths[],
                              const MPI_Aint array_of_displacements[], MPI_Datatype oldtype,
                              MPI_Datatype *newtype);
int PMPI_Type_create_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype,
                             MPI_Datatype *newtype);
int PMPI_Type_create_indexed_block(int count, int blocklength, const int array_of_displacements[],
                                   MPI_Datatype oldtype, MPI_Datatype *newtype);
int PMPI_Type_create_hindexed_block(int count, int blocklength,
                                    const MPI_Aint array_of_displacements[],
                                    MPI_Datatype oldtype, MPI_Datatype *newtype);
int PMPI_Type_create_resized(MPI_Datatype oldtype, MPI_Aint lb, MPI_Aint extent,
                             MPI_Datatype *newtype);
int PMPI_Type_create_struct(int count, const int array_of_blocklengths[],
                            const MPI_Aint array_of_displacements[],
                            const MPI_Datatype array_of_types[], MPI_Datatype *newtype);
int PMPI_Type_create_subarray(int ndims, const int array_of_sizes[],
                              const int array_of_subsizes[], const int array_of_starts[],
                              int order, MPI_Datatype oldtype, MPI_Datatype *newtype);
int PMPI_Type_get_extent(MPI_Datatype datatype, MPI_Aint *lb, MPI_Aint *extent);
int PMPI_Type_get_true_extent(MPI_Datatype datatype, MPI_Aint *true_lb, MPI_Aint *true_extent);
int PMPI_Unpack_external(const char datarep[], const void *inbuf, MPI_Aint insize,
                         MPI_Aint *position, void *outbuf, int outcount, MPI_Datatype datatype)
                         ;
int PMPI_Win_create_errhandler(MPI_Win_errhandler_function *win_errhandler_fn,
                               MPI_Errhandler *errhandler);
int PMPI_Win_get_errhandler(MPI_Win win, MPI_Errhandler *errhandler);
int PMPI_Win_set_errhandler(MPI_Win win, MPI_Errhandler errhandler);




int PMPI_Type_create_f90_integer(int r, MPI_Datatype *newtype);
int PMPI_Type_create_f90_real(int p, int r, MPI_Datatype *newtype);
int PMPI_Type_create_f90_complex(int p, int r, MPI_Datatype *newtype);

int PMPI_Reduce_local(const void *inbuf, void *inoutbuf, int count, MPI_Datatype datatype,
                      MPI_Op op)
                      ;
int PMPI_Op_commutative(MPI_Op op, int *commute);
int PMPI_Reduce_scatter_block(const void *sendbuf, void *recvbuf, int recvcount,
                              MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                             
                              ;
int PMPI_Dist_graph_create_adjacent(MPI_Comm comm_old, int indegree, const int sources[],
                                    const int sourceweights[], int outdegree,
                                    const int destinations[], const int destweights[],
                                    MPI_Info info, int reorder, MPI_Comm *comm_dist_graph);
int PMPI_Dist_graph_create(MPI_Comm comm_old, int n, const int sources[], const int degrees[],
                           const int destinations[], const int weights[], MPI_Info info,
                           int reorder, MPI_Comm *comm_dist_graph);
int PMPI_Dist_graph_neighbors_count(MPI_Comm comm, int *indegree, int *outdegree, int *weighted);
int PMPI_Dist_graph_neighbors(MPI_Comm comm, int maxindegree, int sources[], int sourceweights[],
                              int maxoutdegree, int destinations[], int destweights[]);


int PMPI_Improbe(int source, int tag, MPI_Comm comm, int *flag, MPI_Message *message,
                 MPI_Status *status);
int PMPI_Imrecv(void *buf, int count, MPI_Datatype datatype, MPI_Message *message,
                MPI_Request *request) ;
int PMPI_Mprobe(int source, int tag, MPI_Comm comm, MPI_Message *message, MPI_Status *status);
int PMPI_Mrecv(void *buf, int count, MPI_Datatype datatype, MPI_Message *message,
               MPI_Status *status) ;


int PMPI_Comm_idup(MPI_Comm comm, MPI_Comm *newcomm, MPI_Request *request);
int PMPI_Ibarrier(MPI_Comm comm, MPI_Request *request);
int PMPI_Ibcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm,
                MPI_Request *request) ;
int PMPI_Igather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                 int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm,
                 MPI_Request *request)
                 ;
int PMPI_Igatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                  const int recvcounts[], const int displs[], MPI_Datatype recvtype, int root,
                  MPI_Comm comm, MPI_Request *request)
                  ;
int PMPI_Iscatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                  int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm,
                  MPI_Request *request)
                  ;
int PMPI_Iscatterv(const void *sendbuf, const int sendcounts[], const int displs[],
                   MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                   int root, MPI_Comm comm, MPI_Request *request)
                   ;
int PMPI_Iallgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                    int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                    ;
int PMPI_Iallgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                     const int recvcounts[], const int displs[], MPI_Datatype recvtype,
                     MPI_Comm comm, MPI_Request *request)
                     ;
int PMPI_Ialltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                   int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                   ;
int PMPI_Ialltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[],
                    MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                    const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm,
                    MPI_Request *request)
                    ;
int PMPI_Ialltoallw(const void *sendbuf, const int sendcounts[], const int sdispls[],
                    const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
                    const int rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm,
                    MPI_Request *request);
int PMPI_Ireduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                 MPI_Op op, int root, MPI_Comm comm, MPI_Request *request)
                 ;
int PMPI_Iallreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                    MPI_Op op, MPI_Comm comm, MPI_Request *request)
                    ;
int PMPI_Ireduce_scatter(const void *sendbuf, void *recvbuf, const int recvcounts[],
                         MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request)
                        
                         ;
int PMPI_Ireduce_scatter_block(const void *sendbuf, void *recvbuf, int recvcount,
                               MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                               MPI_Request *request)
                              
                               ;
int PMPI_Iscan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
               MPI_Comm comm, MPI_Request *request)
               ;
int PMPI_Iexscan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                 MPI_Op op, MPI_Comm comm, MPI_Request *request)
                 ;


int PMPI_Ineighbor_allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                             void *recvbuf, int recvcount, MPI_Datatype recvtype,
                             MPI_Comm comm, MPI_Request *request)
                            
                             ;
int PMPI_Ineighbor_allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                              void *recvbuf, const int recvcounts[], const int displs[],
                              MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                             
                              ;
int PMPI_Ineighbor_alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                            void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                            MPI_Request *request)
                           
                            ;
int PMPI_Ineighbor_alltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[],
                             MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                             const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm,
                             MPI_Request *request)
                            
                             ;
int PMPI_Ineighbor_alltoallw(const void *sendbuf, const int sendcounts[],
                             const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
                             void *recvbuf, const int recvcounts[], const MPI_Aint rdispls[],
                             const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Request *request);
int PMPI_Neighbor_allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                            void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                           
                            ;
int PMPI_Neighbor_allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                             void *recvbuf, const int recvcounts[], const int displs[],
                             MPI_Datatype recvtype, MPI_Comm comm)
                            
                             ;
int PMPI_Neighbor_alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                           void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                          
                           ;
int PMPI_Neighbor_alltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[],
                            MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                            const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm)
                           
                            ;
int PMPI_Neighbor_alltoallw(const void *sendbuf, const int sendcounts[], const MPI_Aint sdispls[],
                            const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
                            const MPI_Aint rdispls[], const MPI_Datatype recvtypes[],
                            MPI_Comm comm);


int PMPI_Comm_split_type(MPI_Comm comm, int split_type, int key, MPI_Info info, MPI_Comm *newcomm);


int PMPI_Comm_create_group(MPI_Comm comm, MPI_Group group, int tag, MPI_Comm *newcomm);


int PMPI_Get_elements_x(const MPI_Status *status, MPI_Datatype datatype, MPI_Count *count);
int PMPI_Status_set_elements_x(MPI_Status *status, MPI_Datatype datatype, MPI_Count count);
int PMPI_Type_get_extent_x(MPI_Datatype datatype, MPI_Count *lb, MPI_Count *extent);
int PMPI_Type_get_true_extent_x(MPI_Datatype datatype, MPI_Count *lb, MPI_Count *extent);
int PMPI_Type_size_x(MPI_Datatype datatype, MPI_Count *size);


MPI_Aint PMPI_Aint_add(MPI_Aint base, MPI_Aint disp);
MPI_Aint PMPI_Aint_diff(MPI_Aint addr1, MPI_Aint addr2);





int PMPI_T_init_thread(int required, int *provided);
int PMPI_T_finalize(void);
int PMPI_T_enum_get_info(MPI_T_enum enumtype, int *num, char *name, int *name_len);
int PMPI_T_enum_get_item(MPI_T_enum enumtype, int indx, int *value, char *name, int *name_len);
int PMPI_T_cvar_get_num(int *num_cvar);
int PMPI_T_cvar_get_info(int cvar_index, char *name, int *name_len, int *verbosity,
                         MPI_Datatype *datatype, MPI_T_enum *enumtype, char *desc, int *desc_len,
                         int *binding, int *scope);
int PMPI_T_cvar_handle_alloc(int cvar_index, void *obj_handle, MPI_T_cvar_handle *handle,
                             int *count);
int PMPI_T_cvar_handle_free(MPI_T_cvar_handle *handle);
int PMPI_T_cvar_read(MPI_T_cvar_handle handle, void *buf);
int PMPI_T_cvar_write(MPI_T_cvar_handle handle, const void *buf);
int PMPI_T_pvar_get_num(int *num_pvar);
int PMPI_T_pvar_get_info(int pvar_index, char *name, int *name_len, int *verbosity, int *var_class,
                         MPI_Datatype *datatype, MPI_T_enum *enumtype, char *desc, int *desc_len,
                         int *binding, int *readonly, int *continuous, int *atomic);
int PMPI_T_pvar_session_create(MPI_T_pvar_session *session);
int PMPI_T_pvar_session_free(MPI_T_pvar_session *session);
int PMPI_T_pvar_handle_alloc(MPI_T_pvar_session session, int pvar_index, void *obj_handle,
                             MPI_T_pvar_handle *handle, int *count);
int PMPI_T_pvar_handle_free(MPI_T_pvar_session session, MPI_T_pvar_handle *handle);
int PMPI_T_pvar_start(MPI_T_pvar_session session, MPI_T_pvar_handle handle);
int PMPI_T_pvar_stop(MPI_T_pvar_session session, MPI_T_pvar_handle handle);
int PMPI_T_pvar_read(MPI_T_pvar_session session, MPI_T_pvar_handle handle, void *buf);
int PMPI_T_pvar_write(MPI_T_pvar_session session, MPI_T_pvar_handle handle, const void *buf);
int PMPI_T_pvar_reset(MPI_T_pvar_session session, MPI_T_pvar_handle handle);
int PMPI_T_pvar_readreset(MPI_T_pvar_session session, MPI_T_pvar_handle handle, void *buf);
int PMPI_T_category_get_num(int *num_cat);
int PMPI_T_category_get_info(int cat_index, char *name, int *name_len, char *desc, int *desc_len,
                             int *num_cvars, int *num_pvars, int *num_categories);
int PMPI_T_category_get_cvars(int cat_index, int len, int indices[]);
int PMPI_T_category_get_pvars(int cat_index, int len, int indices[]);
int PMPI_T_category_get_categories(int cat_index, int len, int indices[]);
int PMPI_T_category_changed(int *stamp);
int PMPI_T_cvar_get_index(const char *name, int *cvar_index);
int PMPI_T_pvar_get_index(const char *name, int var_class, int *pvar_index);
int PMPI_T_category_get_index(const char *name, int *cat_index);





int PMPIX_Comm_failure_ack(MPI_Comm comm);
int PMPIX_Comm_failure_get_acked(MPI_Comm comm, MPI_Group *failedgrp);
int PMPIX_Comm_revoke(MPI_Comm comm);
int PMPIX_Comm_shrink(MPI_Comm comm, MPI_Comm *newcomm);
int PMPIX_Comm_agree(MPI_Comm comm, int *flag);
# 2315 "/sw/eb/sw/impi/2018.5.288-iccifort-2019.5.281/include64/mpi.h"
# 1 "/sw/eb/sw/impi/2018.5.288-iccifort-2019.5.281/include64/mpio.h" 1
# 77 "/sw/eb/sw/impi/2018.5.288-iccifort-2019.5.281/include64/mpio.h"
# 1 "/sw/eb/sw/impi/2018.5.288-iccifort-2019.5.281/include64/mpi.h" 1
# 78 "/sw/eb/sw/impi/2018.5.288-iccifort-2019.5.281/include64/mpio.h" 2
# 187 "/sw/eb/sw/impi/2018.5.288-iccifort-2019.5.281/include64/mpio.h"
int MPI_File_open(MPI_Comm comm, const char *filename, int amode, MPI_Info info, MPI_File *fh);
int MPI_File_close(MPI_File *fh);
int MPI_File_delete(const char *filename, MPI_Info info);
int MPI_File_set_size(MPI_File fh, MPI_Offset size);
int MPI_File_preallocate(MPI_File fh, MPI_Offset size);
int MPI_File_get_size(MPI_File fh, MPI_Offset *size);
int MPI_File_get_group(MPI_File fh, MPI_Group *group);
int MPI_File_get_amode(MPI_File fh, int *amode);
int MPI_File_set_info(MPI_File fh, MPI_Info info);
int MPI_File_get_info(MPI_File fh, MPI_Info *info_used);


int MPI_File_set_view(MPI_File fh, MPI_Offset disp, MPI_Datatype etype, MPI_Datatype filetype,
                      const char *datarep, MPI_Info info);
int MPI_File_get_view(MPI_File fh, MPI_Offset *disp, MPI_Datatype *etype, MPI_Datatype *filetype,
                      char *datarep);


int MPI_File_read_at(MPI_File fh, MPI_Offset offset, void *buf, int count, MPI_Datatype datatype,
                     MPI_Status *status) ;
int MPI_File_read_at_all(MPI_File fh, MPI_Offset offset, void * buf, int count,
                         MPI_Datatype datatype, MPI_Status *status)
    ;
int MPI_File_write_at(MPI_File fh, MPI_Offset offset, const void * buf, int count,
                      MPI_Datatype datatype, MPI_Status *status)
    ;
int MPI_File_write_at_all(MPI_File fh, MPI_Offset offset, const void *buf, int count,
                          MPI_Datatype datatype, MPI_Status *status)
    ;




int MPI_File_iread_at(MPI_File fh, MPI_Offset offset, void *buf, int count, MPI_Datatype datatype,
                      MPI_Request *request) ;
int MPI_File_iwrite_at(MPI_File fh, MPI_Offset offset, const void *buf, int count,
                       MPI_Datatype datatype, MPI_Request *request)
    ;


int MPI_File_read(MPI_File fh, void *buf, int count, MPI_Datatype datatype, MPI_Status *status)
    ;
int MPI_File_read_all(MPI_File fh, void *buf, int count, MPI_Datatype datatype, MPI_Status *status)
    ;
int MPI_File_write(MPI_File fh, const void *buf, int count, MPI_Datatype datatype,
                   MPI_Status *status) ;
int MPI_File_write_all(MPI_File fh, const void *buf, int count, MPI_Datatype datatype,
                       MPI_Status *status) ;





int MPI_File_iread(MPI_File fh, void *buf, int count, MPI_Datatype datatype, MPI_Request *request)
    ;
int MPI_File_iwrite(MPI_File fh, const void *buf, int count, MPI_Datatype datatype,
                    MPI_Request *request) ;

int MPI_File_seek(MPI_File fh, MPI_Offset offset, int whence);
int MPI_File_get_position(MPI_File fh, MPI_Offset *offset);
int MPI_File_get_byte_offset(MPI_File fh, MPI_Offset offset, MPI_Offset *disp);


int MPI_File_read_shared(MPI_File fh, void *buf, int count, MPI_Datatype datatype,
                         MPI_Status *status) ;
int MPI_File_write_shared(MPI_File fh, const void *buf, int count, MPI_Datatype datatype,
                          MPI_Status *status) ;
int MPI_File_iread_shared(MPI_File fh, void *buf, int count, MPI_Datatype datatype,
                          MPI_Request *request) ;
int MPI_File_iwrite_shared(MPI_File fh, const void *buf, int count, MPI_Datatype datatype,
                           MPI_Request *request) ;
int MPI_File_read_ordered(MPI_File fh, void *buf, int count, MPI_Datatype datatype,
                          MPI_Status *status) ;
int MPI_File_write_ordered(MPI_File fh, const void *buf, int count, MPI_Datatype datatype,
                           MPI_Status *status) ;
int MPI_File_seek_shared(MPI_File fh, MPI_Offset offset, int whence);
int MPI_File_get_position_shared(MPI_File fh, MPI_Offset *offset);


int MPI_File_read_at_all_begin(MPI_File fh, MPI_Offset offset, void *buf, int count,
                               MPI_Datatype datatype) ;
int MPI_File_read_at_all_end(MPI_File fh, void *buf, MPI_Status *status);
int MPI_File_write_at_all_begin(MPI_File fh, MPI_Offset offset, const void *buf, int count,
                                MPI_Datatype datatype) ;
int MPI_File_write_at_all_end(MPI_File fh, const void *buf, MPI_Status *status);
int MPI_File_read_all_begin(MPI_File fh, void *buf, int count, MPI_Datatype datatype)
    ;
int MPI_File_read_all_end(MPI_File fh, void *buf, MPI_Status *status);
int MPI_File_write_all_begin(MPI_File fh, const void *buf, int count, MPI_Datatype datatype)
    ;
int MPI_File_write_all_end(MPI_File fh, const void *buf, MPI_Status *status);
int MPI_File_read_ordered_begin(MPI_File fh, void *buf, int count, MPI_Datatype datatype)
    ;
int MPI_File_read_ordered_end(MPI_File fh, void *buf, MPI_Status *status);
int MPI_File_write_ordered_begin(MPI_File fh, const void *buf, int count, MPI_Datatype datatype)
    ;
int MPI_File_write_ordered_end(MPI_File fh, const void *buf, MPI_Status *status);


int MPI_File_get_type_extent(MPI_File fh, MPI_Datatype datatype, MPI_Aint *extent);


int MPI_Register_datarep(const char *datarep, MPI_Datarep_conversion_function *read_conversion_fn,
    MPI_Datarep_conversion_function *write_conversion_fn,
    MPI_Datarep_extent_function *dtype_file_extent_fn, void *extra_state);


int MPI_File_set_atomicity(MPI_File fh, int flag);
int MPI_File_get_atomicity(MPI_File fh, int *flag);
int MPI_File_sync(MPI_File fh);
# 306 "/sw/eb/sw/impi/2018.5.288-iccifort-2019.5.281/include64/mpio.h"
int MPI_File_iread_at_all(MPI_File fh, MPI_Offset offset, void *buf, int count,
                           MPI_Datatype datatype, MPI_Request *request)
    ;
int MPI_File_iwrite_at_all(MPI_File fh, MPI_Offset offset, const void *buf, int count,
                            MPI_Datatype datatype, MPI_Request *request)
    ;
int MPI_File_iread_all(MPI_File fh, void *buf, int count, MPI_Datatype datatype,
                        MPI_Request *request)
    ;
int MPI_File_iwrite_all(MPI_File fh, const void *buf, int count, MPI_Datatype datatype,
                         MPI_Request *request)
    ;
# 346 "/sw/eb/sw/impi/2018.5.288-iccifort-2019.5.281/include64/mpio.h"
MPI_File MPI_File_f2c(MPI_Fint file);
MPI_Fint MPI_File_c2f(MPI_File file);
# 407 "/sw/eb/sw/impi/2018.5.288-iccifort-2019.5.281/include64/mpio.h"
int PMPI_File_open(MPI_Comm, const char *, int, MPI_Info, MPI_File *);
int PMPI_File_close(MPI_File *);
int PMPI_File_delete(const char *, MPI_Info);
int PMPI_File_set_size(MPI_File, MPI_Offset);
int PMPI_File_preallocate(MPI_File, MPI_Offset);
int PMPI_File_get_size(MPI_File, MPI_Offset *);
int PMPI_File_get_group(MPI_File, MPI_Group *);
int PMPI_File_get_amode(MPI_File, int *);
int PMPI_File_set_info(MPI_File, MPI_Info);
int PMPI_File_get_info(MPI_File, MPI_Info *);


int PMPI_File_set_view(MPI_File, MPI_Offset,
    MPI_Datatype, MPI_Datatype, const char *, MPI_Info);
int PMPI_File_get_view(MPI_File, MPI_Offset *,
      MPI_Datatype *, MPI_Datatype *, char *);


int PMPI_File_read_at(MPI_File, MPI_Offset, void *,
       int, MPI_Datatype, MPI_Status *)
              ;
int PMPI_File_read_at_all(MPI_File, MPI_Offset, void *,
       int, MPI_Datatype, MPI_Status *)
              ;
int PMPI_File_write_at(MPI_File, MPI_Offset, const void *,
       int, MPI_Datatype, MPI_Status *)
              ;
int PMPI_File_write_at_all(MPI_File, MPI_Offset, const void *,
       int, MPI_Datatype, MPI_Status *)
              ;





int PMPI_File_iread_at(MPI_File, MPI_Offset, void *,
       int, MPI_Datatype, MPI_Request *)
              ;
int PMPI_File_iwrite_at(MPI_File, MPI_Offset, const void *,
       int, MPI_Datatype, MPI_Request *)
              ;


int PMPI_File_read(MPI_File, void *, int, MPI_Datatype, MPI_Status *)
                   ;
int PMPI_File_read_all(MPI_File, void *, int, MPI_Datatype, MPI_Status *)
                       ;
int PMPI_File_write(MPI_File, const void *, int, MPI_Datatype, MPI_Status *)
                    ;
int PMPI_File_write_all(MPI_File, const void *, int, MPI_Datatype, MPI_Status *)
                        ;





int PMPI_File_iread(MPI_File, void *, int, MPI_Datatype, MPI_Request *)
                    ;
int PMPI_File_iwrite(MPI_File, const void *, int, MPI_Datatype, MPI_Request *)
                     ;

int PMPI_File_seek(MPI_File, MPI_Offset, int);
int PMPI_File_get_position(MPI_File, MPI_Offset *);
int PMPI_File_get_byte_offset(MPI_File, MPI_Offset, MPI_Offset *);


int PMPI_File_read_shared(MPI_File, void *, int, MPI_Datatype, MPI_Status *)
                          ;
int PMPI_File_write_shared(MPI_File, const void *, int, MPI_Datatype, MPI_Status *)
                           ;
int PMPI_File_iread_shared(MPI_File, void *, int,
      MPI_Datatype, MPI_Request *)
                           ;
int PMPI_File_iwrite_shared(MPI_File, const void *, int,
       MPI_Datatype, MPI_Request *)
                            ;
int PMPI_File_read_ordered(MPI_File, void *, int, MPI_Datatype, MPI_Status *)
                           ;
int PMPI_File_write_ordered(MPI_File, const void *, int, MPI_Datatype, MPI_Status *)
                            ;
int PMPI_File_seek_shared(MPI_File, MPI_Offset, int);
int PMPI_File_get_position_shared(MPI_File, MPI_Offset *);


int PMPI_File_read_at_all_begin(MPI_File, MPI_Offset, void *,
                               int, MPI_Datatype)
                               ;
int PMPI_File_read_at_all_end(MPI_File, void *, MPI_Status *);
int PMPI_File_write_at_all_begin(MPI_File, MPI_Offset, const void *,
                                 int, MPI_Datatype)
                                 ;
int PMPI_File_write_at_all_end(MPI_File, const void *, MPI_Status *);
int PMPI_File_read_all_begin(MPI_File, void *, int, MPI_Datatype)
                             ;
int PMPI_File_read_all_end(MPI_File, void *, MPI_Status *);
int PMPI_File_write_all_begin(MPI_File, const void *, int, MPI_Datatype)
                              ;
int PMPI_File_write_all_end(MPI_File, const void *, MPI_Status *);
int PMPI_File_read_ordered_begin(MPI_File, void *, int, MPI_Datatype)
                                 ;
int PMPI_File_read_ordered_end(MPI_File, void *, MPI_Status *);
int PMPI_File_write_ordered_begin(MPI_File, const void *, int, MPI_Datatype)
                                  ;
int PMPI_File_write_ordered_end(MPI_File, const void *, MPI_Status *);


int PMPI_File_get_type_extent(MPI_File, MPI_Datatype, MPI_Aint *);


int PMPI_Register_datarep(const char *,
    MPI_Datarep_conversion_function *,
    MPI_Datarep_conversion_function *,
    MPI_Datarep_extent_function *,
    void *);


int PMPI_File_set_atomicity(MPI_File, int);
int PMPI_File_get_atomicity(MPI_File, int *);
int PMPI_File_sync(MPI_File);
# 535 "/sw/eb/sw/impi/2018.5.288-iccifort-2019.5.281/include64/mpio.h"
int PMPI_File_iread_at_all(MPI_File fh, MPI_Offset offset, void *buf, int count,
                            MPI_Datatype datatype, MPI_Request *request)
    ;
int PMPI_File_iwrite_at_all(MPI_File fh, MPI_Offset offset, const void *buf, int count,
                             MPI_Datatype datatype, MPI_Request *request)
    ;
int PMPI_File_iread_all(MPI_File fh, void *buf, int count, MPI_Datatype datatype,
                         MPI_Request *request)
    ;
int PMPI_File_iwrite_all(MPI_File fh, const void *buf, int count, MPI_Datatype datatype,
                          MPI_Request *request)
    ;
# 559 "/sw/eb/sw/impi/2018.5.288-iccifort-2019.5.281/include64/mpio.h"
MPI_File PMPI_File_f2c(MPI_Fint);
MPI_Fint PMPI_File_c2f(MPI_File);
# 2316 "/sw/eb/sw/impi/2018.5.288-iccifort-2019.5.281/include64/mpi.h" 2
# 2337 "/sw/eb/sw/impi/2018.5.288-iccifort-2019.5.281/include64/mpi.h"
typedef int MPIX_Grequest_class;
int MPIX_Grequest_class_create(MPI_Grequest_query_function *query_fn,
                               MPI_Grequest_free_function *free_fn,
                               MPI_Grequest_cancel_function *cancel_fn,
                               MPIX_Grequest_poll_function *poll_fn,
                               MPIX_Grequest_wait_function *wait_fn,
                               MPIX_Grequest_class *greq_class);
int MPIX_Grequest_class_allocate(MPIX_Grequest_class greq_class, void *extra_state,
                                 MPI_Request *request);
int MPIX_Grequest_start(MPI_Grequest_query_function *query_fn,
                        MPI_Grequest_free_function *free_fn,
                        MPI_Grequest_cancel_function *cancel_fn,
                        MPIX_Grequest_poll_function *poll_fn,
                        MPIX_Grequest_wait_function *wait_fn, void *extra_state,
                        MPI_Request *request);


struct mpixi_mutex_s;
typedef struct mpixi_mutex_s * MPIX_Mutex;
int MPIX_Mutex_create(int count, MPI_Comm comm, MPIX_Mutex *hdl);
int MPIX_Mutex_free(MPIX_Mutex *hdl);
int MPIX_Mutex_lock(MPIX_Mutex hdl, int mutex, int proc);
int MPIX_Mutex_unlock(MPIX_Mutex hdl, int mutex, int proc);




int PMPIX_Grequest_class_create(MPI_Grequest_query_function *query_fn,
                                MPI_Grequest_free_function *free_fn,
                                MPI_Grequest_cancel_function *cancel_fn,
                                MPIX_Grequest_poll_function *poll_fn,
                                MPIX_Grequest_wait_function *wait_fn,
                                MPIX_Grequest_class *greq_class);
int PMPIX_Grequest_class_allocate(MPIX_Grequest_class greq_class, void *extra_state,
                                  MPI_Request *request);
int PMPIX_Grequest_start(MPI_Grequest_query_function *query_fn,
                         MPI_Grequest_free_function *free_fn,
                         MPI_Grequest_cancel_function *cancel_fn,
                         MPIX_Grequest_poll_function *poll_fn,
                         MPIX_Grequest_wait_function *wait_fn, void *extra_state,
                         MPI_Request *request);


int PMPIX_Mutex_create(int count, MPI_Comm comm, MPIX_Mutex *hdl);
int PMPIX_Mutex_free(MPIX_Mutex *hdl);
int PMPIX_Mutex_lock(MPIX_Mutex hdl, int mutex, int proc);
int PMPIX_Mutex_unlock(MPIX_Mutex hdl, int mutex, int proc);
# 109 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h" 2
# 138 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
# 1 "/usr/include/stdio.h" 1 3 4
# 29 "/usr/include/stdio.h" 3 4




# 1 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/stddef.h" 1 3 4
# 216 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/stddef.h" 3 4

# 216 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/stddef.h" 3 4
typedef long unsigned int size_t;
# 34 "/usr/include/stdio.h" 2 3 4

# 1 "/usr/include/bits/types.h" 1 3 4
# 27 "/usr/include/bits/types.h" 3 4
# 1 "/usr/include/bits/wordsize.h" 1 3 4
# 28 "/usr/include/bits/types.h" 2 3 4


typedef unsigned char __u_char;
typedef unsigned short int __u_short;
typedef unsigned int __u_int;
typedef unsigned long int __u_long;


typedef signed char __int8_t;
typedef unsigned char __uint8_t;
typedef signed short int __int16_t;
typedef unsigned short int __uint16_t;
typedef signed int __int32_t;
typedef unsigned int __uint32_t;

typedef signed long int __int64_t;
typedef unsigned long int __uint64_t;







typedef long int __quad_t;
typedef unsigned long int __u_quad_t;
# 130 "/usr/include/bits/types.h" 3 4
# 1 "/usr/include/bits/typesizes.h" 1 3 4
# 131 "/usr/include/bits/types.h" 2 3 4


typedef unsigned long int __dev_t;
typedef unsigned int __uid_t;
typedef unsigned int __gid_t;
typedef unsigned long int __ino_t;
typedef unsigned long int __ino64_t;
typedef unsigned int __mode_t;
typedef unsigned long int __nlink_t;
typedef long int __off_t;
typedef long int __off64_t;
typedef int __pid_t;
typedef struct { int __val[2]; } __fsid_t;
typedef long int __clock_t;
typedef unsigned long int __rlim_t;
typedef unsigned long int __rlim64_t;
typedef unsigned int __id_t;
typedef long int __time_t;
typedef unsigned int __useconds_t;
typedef long int __suseconds_t;

typedef int __daddr_t;
typedef int __key_t;


typedef int __clockid_t;


typedef void * __timer_t;


typedef long int __blksize_t;




typedef long int __blkcnt_t;
typedef long int __blkcnt64_t;


typedef unsigned long int __fsblkcnt_t;
typedef unsigned long int __fsblkcnt64_t;


typedef unsigned long int __fsfilcnt_t;
typedef unsigned long int __fsfilcnt64_t;


typedef long int __fsword_t;

typedef long int __ssize_t;


typedef long int __syscall_slong_t;

typedef unsigned long int __syscall_ulong_t;



typedef __off64_t __loff_t;
typedef __quad_t *__qaddr_t;
typedef char *__caddr_t;


typedef long int __intptr_t;


typedef unsigned int __socklen_t;
# 36 "/usr/include/stdio.h" 2 3 4
# 44 "/usr/include/stdio.h" 3 4
struct _IO_FILE;



typedef struct _IO_FILE FILE;





# 64 "/usr/include/stdio.h" 3 4
typedef struct _IO_FILE __FILE;
# 74 "/usr/include/stdio.h" 3 4
# 1 "/usr/include/libio.h" 1 3 4
# 32 "/usr/include/libio.h" 3 4
# 1 "/usr/include/_G_config.h" 1 3 4
# 15 "/usr/include/_G_config.h" 3 4
# 1 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/stddef.h" 1 3 4
# 16 "/usr/include/_G_config.h" 2 3 4




# 1 "/usr/include/wchar.h" 1 3 4
# 82 "/usr/include/wchar.h" 3 4
typedef struct
{
  int __count;
  union
  {

    unsigned int __wch;



    char __wchb[4];
  } __value;
} __mbstate_t;
# 21 "/usr/include/_G_config.h" 2 3 4
typedef struct
{
  __off_t __pos;
  __mbstate_t __state;
} _G_fpos_t;
typedef struct
{
  __off64_t __pos;
  __mbstate_t __state;
} _G_fpos64_t;
# 33 "/usr/include/libio.h" 2 3 4
# 50 "/usr/include/libio.h" 3 4
# 1 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/stdarg.h" 1 3 4
# 40 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/stdarg.h" 3 4
typedef __builtin_va_list __gnuc_va_list;
# 51 "/usr/include/libio.h" 2 3 4
# 145 "/usr/include/libio.h" 3 4
struct _IO_jump_t; struct _IO_FILE;
# 155 "/usr/include/libio.h" 3 4
typedef void _IO_lock_t;





struct _IO_marker {
  struct _IO_marker *_next;
  struct _IO_FILE *_sbuf;



  int _pos;
# 178 "/usr/include/libio.h" 3 4
};


enum __codecvt_result
{
  __codecvt_ok,
  __codecvt_partial,
  __codecvt_error,
  __codecvt_noconv
};
# 246 "/usr/include/libio.h" 3 4
struct _IO_FILE {
  int _flags;




  char* _IO_read_ptr;
  char* _IO_read_end;
  char* _IO_read_base;
  char* _IO_write_base;
  char* _IO_write_ptr;
  char* _IO_write_end;
  char* _IO_buf_base;
  char* _IO_buf_end;

  char *_IO_save_base;
  char *_IO_backup_base;
  char *_IO_save_end;

  struct _IO_marker *_markers;

  struct _IO_FILE *_chain;

  int _fileno;



  int _flags2;

  __off_t _old_offset;



  unsigned short _cur_column;
  signed char _vtable_offset;
  char _shortbuf[1];



  _IO_lock_t *_lock;
# 294 "/usr/include/libio.h" 3 4
  __off64_t _offset;
# 303 "/usr/include/libio.h" 3 4
  void *__pad1;
  void *__pad2;
  void *__pad3;
  void *__pad4;
  size_t __pad5;

  int _mode;

  char _unused2[15 * sizeof (int) - 4 * sizeof (void *) - sizeof (size_t)];

};


typedef struct _IO_FILE _IO_FILE;


struct _IO_FILE_plus;

extern struct _IO_FILE_plus _IO_2_1_stdin_;
extern struct _IO_FILE_plus _IO_2_1_stdout_;
extern struct _IO_FILE_plus _IO_2_1_stderr_;
# 339 "/usr/include/libio.h" 3 4
typedef __ssize_t __io_read_fn (void *__cookie, char *__buf, size_t __nbytes);







typedef __ssize_t __io_write_fn (void *__cookie, const char *__buf,
     size_t __n);







typedef int __io_seek_fn (void *__cookie, __off64_t *__pos, int __w);


typedef int __io_close_fn (void *__cookie);
# 391 "/usr/include/libio.h" 3 4
extern int __underflow (_IO_FILE *);
extern int __uflow (_IO_FILE *);
extern int __overflow (_IO_FILE *, int);
# 435 "/usr/include/libio.h" 3 4
extern int _IO_getc (_IO_FILE *__fp);
extern int _IO_putc (int __c, _IO_FILE *__fp);
extern int _IO_feof (_IO_FILE *__fp) __attribute__ ((__nothrow__ , __leaf__));
extern int _IO_ferror (_IO_FILE *__fp) __attribute__ ((__nothrow__ , __leaf__));

extern int _IO_peekc_locked (_IO_FILE *__fp);





extern void _IO_flockfile (_IO_FILE *) __attribute__ ((__nothrow__ , __leaf__));
extern void _IO_funlockfile (_IO_FILE *) __attribute__ ((__nothrow__ , __leaf__));
extern int _IO_ftrylockfile (_IO_FILE *) __attribute__ ((__nothrow__ , __leaf__));
# 465 "/usr/include/libio.h" 3 4
extern int _IO_vfscanf (_IO_FILE * __restrict, const char * __restrict,
   __gnuc_va_list, int *__restrict);
extern int _IO_vfprintf (_IO_FILE *__restrict, const char *__restrict,
    __gnuc_va_list);
extern __ssize_t _IO_padn (_IO_FILE *, int, __ssize_t);
extern size_t _IO_sgetn (_IO_FILE *, void *, size_t);

extern __off64_t _IO_seekoff (_IO_FILE *, __off64_t, int, int);
extern __off64_t _IO_seekpos (_IO_FILE *, __off64_t, int);

extern void _IO_free_backup_area (_IO_FILE *) __attribute__ ((__nothrow__ , __leaf__));
# 75 "/usr/include/stdio.h" 2 3 4




typedef __gnuc_va_list va_list;
# 90 "/usr/include/stdio.h" 3 4
typedef __off_t off_t;
# 102 "/usr/include/stdio.h" 3 4
typedef __ssize_t ssize_t;







typedef _G_fpos_t fpos_t;




# 164 "/usr/include/stdio.h" 3 4
# 1 "/usr/include/bits/stdio_lim.h" 1 3 4
# 165 "/usr/include/stdio.h" 2 3 4



extern struct _IO_FILE *stdin;
extern struct _IO_FILE *stdout;
extern struct _IO_FILE *stderr;







extern int remove (const char *__filename) __attribute__ ((__nothrow__ , __leaf__));

extern int rename (const char *__old, const char *__new) __attribute__ ((__nothrow__ , __leaf__));




extern int renameat (int __oldfd, const char *__old, int __newfd,
       const char *__new) __attribute__ ((__nothrow__ , __leaf__));








extern FILE *tmpfile (void) ;
# 209 "/usr/include/stdio.h" 3 4
extern char *tmpnam (char *__s) __attribute__ ((__nothrow__ , __leaf__)) ;





extern char *tmpnam_r (char *__s) __attribute__ ((__nothrow__ , __leaf__)) ;
# 227 "/usr/include/stdio.h" 3 4
extern char *tempnam (const char *__dir, const char *__pfx)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__malloc__)) ;








extern int fclose (FILE *__stream);




extern int fflush (FILE *__stream);

# 252 "/usr/include/stdio.h" 3 4
extern int fflush_unlocked (FILE *__stream);
# 266 "/usr/include/stdio.h" 3 4






extern FILE *fopen (const char *__restrict __filename,
      const char *__restrict __modes) ;




extern FILE *freopen (const char *__restrict __filename,
        const char *__restrict __modes,
        FILE *__restrict __stream) ;
# 295 "/usr/include/stdio.h" 3 4

# 306 "/usr/include/stdio.h" 3 4
extern FILE *fdopen (int __fd, const char *__modes) __attribute__ ((__nothrow__ , __leaf__)) ;
# 319 "/usr/include/stdio.h" 3 4
extern FILE *fmemopen (void *__s, size_t __len, const char *__modes)
  __attribute__ ((__nothrow__ , __leaf__)) ;




extern FILE *open_memstream (char **__bufloc, size_t *__sizeloc) __attribute__ ((__nothrow__ , __leaf__)) ;






extern void setbuf (FILE *__restrict __stream, char *__restrict __buf) __attribute__ ((__nothrow__ , __leaf__));



extern int setvbuf (FILE *__restrict __stream, char *__restrict __buf,
      int __modes, size_t __n) __attribute__ ((__nothrow__ , __leaf__));





extern void setbuffer (FILE *__restrict __stream, char *__restrict __buf,
         size_t __size) __attribute__ ((__nothrow__ , __leaf__));


extern void setlinebuf (FILE *__stream) __attribute__ ((__nothrow__ , __leaf__));








extern int fprintf (FILE *__restrict __stream,
      const char *__restrict __format, ...);




extern int printf (const char *__restrict __format, ...);

extern int sprintf (char *__restrict __s,
      const char *__restrict __format, ...) __attribute__ ((__nothrow__));





extern int vfprintf (FILE *__restrict __s, const char *__restrict __format,
       __gnuc_va_list __arg);




extern int vprintf (const char *__restrict __format, __gnuc_va_list __arg);

extern int vsprintf (char *__restrict __s, const char *__restrict __format,
       __gnuc_va_list __arg) __attribute__ ((__nothrow__));





extern int snprintf (char *__restrict __s, size_t __maxlen,
       const char *__restrict __format, ...)
     __attribute__ ((__nothrow__)) __attribute__ ((__format__ (__printf__, 3, 4)));

extern int vsnprintf (char *__restrict __s, size_t __maxlen,
        const char *__restrict __format, __gnuc_va_list __arg)
     __attribute__ ((__nothrow__)) __attribute__ ((__format__ (__printf__, 3, 0)));

# 412 "/usr/include/stdio.h" 3 4
extern int vdprintf (int __fd, const char *__restrict __fmt,
       __gnuc_va_list __arg)
     __attribute__ ((__format__ (__printf__, 2, 0)));
extern int dprintf (int __fd, const char *__restrict __fmt, ...)
     __attribute__ ((__format__ (__printf__, 2, 3)));








extern int fscanf (FILE *__restrict __stream,
     const char *__restrict __format, ...) ;




extern int scanf (const char *__restrict __format, ...) ;

extern int sscanf (const char *__restrict __s,
     const char *__restrict __format, ...) __attribute__ ((__nothrow__ , __leaf__));
# 443 "/usr/include/stdio.h" 3 4
extern int fscanf (FILE *__restrict __stream, const char *__restrict __format, ...) __asm__ ("" "__isoc99_fscanf")

                               ;
extern int scanf (const char *__restrict __format, ...) __asm__ ("" "__isoc99_scanf")
                              ;
extern int sscanf (const char *__restrict __s, const char *__restrict __format, ...) __asm__ ("" "__isoc99_sscanf") __attribute__ ((__nothrow__ , __leaf__))

                      ;
# 463 "/usr/include/stdio.h" 3 4








extern int vfscanf (FILE *__restrict __s, const char *__restrict __format,
      __gnuc_va_list __arg)
     __attribute__ ((__format__ (__scanf__, 2, 0))) ;





extern int vscanf (const char *__restrict __format, __gnuc_va_list __arg)
     __attribute__ ((__format__ (__scanf__, 1, 0))) ;


extern int vsscanf (const char *__restrict __s,
      const char *__restrict __format, __gnuc_va_list __arg)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__format__ (__scanf__, 2, 0)));
# 494 "/usr/include/stdio.h" 3 4
extern int vfscanf (FILE *__restrict __s, const char *__restrict __format, __gnuc_va_list __arg) __asm__ ("" "__isoc99_vfscanf")



     __attribute__ ((__format__ (__scanf__, 2, 0))) ;
extern int vscanf (const char *__restrict __format, __gnuc_va_list __arg) __asm__ ("" "__isoc99_vscanf")

     __attribute__ ((__format__ (__scanf__, 1, 0))) ;
extern int vsscanf (const char *__restrict __s, const char *__restrict __format, __gnuc_va_list __arg) __asm__ ("" "__isoc99_vsscanf") __attribute__ ((__nothrow__ , __leaf__))



     __attribute__ ((__format__ (__scanf__, 2, 0)));
# 522 "/usr/include/stdio.h" 3 4









extern int fgetc (FILE *__stream);
extern int getc (FILE *__stream);





extern int getchar (void);

# 550 "/usr/include/stdio.h" 3 4
extern int getc_unlocked (FILE *__stream);
extern int getchar_unlocked (void);
# 561 "/usr/include/stdio.h" 3 4
extern int fgetc_unlocked (FILE *__stream);











extern int fputc (int __c, FILE *__stream);
extern int putc (int __c, FILE *__stream);





extern int putchar (int __c);

# 594 "/usr/include/stdio.h" 3 4
extern int fputc_unlocked (int __c, FILE *__stream);







extern int putc_unlocked (int __c, FILE *__stream);
extern int putchar_unlocked (int __c);






extern int getw (FILE *__stream);


extern int putw (int __w, FILE *__stream);








extern char *fgets (char *__restrict __s, int __n, FILE *__restrict __stream)
     ;
# 640 "/usr/include/stdio.h" 3 4

# 665 "/usr/include/stdio.h" 3 4
extern __ssize_t __getdelim (char **__restrict __lineptr,
          size_t *__restrict __n, int __delimiter,
          FILE *__restrict __stream) ;
extern __ssize_t getdelim (char **__restrict __lineptr,
        size_t *__restrict __n, int __delimiter,
        FILE *__restrict __stream) ;







extern __ssize_t getline (char **__restrict __lineptr,
       size_t *__restrict __n,
       FILE *__restrict __stream) ;








extern int fputs (const char *__restrict __s, FILE *__restrict __stream);





extern int puts (const char *__s);






extern int ungetc (int __c, FILE *__stream);






extern size_t fread (void *__restrict __ptr, size_t __size,
       size_t __n, FILE *__restrict __stream) ;




extern size_t fwrite (const void *__restrict __ptr, size_t __size,
        size_t __n, FILE *__restrict __s);

# 737 "/usr/include/stdio.h" 3 4
extern size_t fread_unlocked (void *__restrict __ptr, size_t __size,
         size_t __n, FILE *__restrict __stream) ;
extern size_t fwrite_unlocked (const void *__restrict __ptr, size_t __size,
          size_t __n, FILE *__restrict __stream);








extern int fseek (FILE *__stream, long int __off, int __whence);




extern long int ftell (FILE *__stream) ;




extern void rewind (FILE *__stream);

# 773 "/usr/include/stdio.h" 3 4
extern int fseeko (FILE *__stream, __off_t __off, int __whence);




extern __off_t ftello (FILE *__stream) ;
# 792 "/usr/include/stdio.h" 3 4






extern int fgetpos (FILE *__restrict __stream, fpos_t *__restrict __pos);




extern int fsetpos (FILE *__stream, const fpos_t *__pos);
# 815 "/usr/include/stdio.h" 3 4

# 824 "/usr/include/stdio.h" 3 4


extern void clearerr (FILE *__stream) __attribute__ ((__nothrow__ , __leaf__));

extern int feof (FILE *__stream) __attribute__ ((__nothrow__ , __leaf__)) ;

extern int ferror (FILE *__stream) __attribute__ ((__nothrow__ , __leaf__)) ;




extern void clearerr_unlocked (FILE *__stream) __attribute__ ((__nothrow__ , __leaf__));
extern int feof_unlocked (FILE *__stream) __attribute__ ((__nothrow__ , __leaf__)) ;
extern int ferror_unlocked (FILE *__stream) __attribute__ ((__nothrow__ , __leaf__)) ;








extern void perror (const char *__s);






# 1 "/usr/include/bits/sys_errlist.h" 1 3 4
# 26 "/usr/include/bits/sys_errlist.h" 3 4
extern int sys_nerr;
extern const char *const sys_errlist[];
# 854 "/usr/include/stdio.h" 2 3 4




extern int fileno (FILE *__stream) __attribute__ ((__nothrow__ , __leaf__)) ;




extern int fileno_unlocked (FILE *__stream) __attribute__ ((__nothrow__ , __leaf__)) ;
# 873 "/usr/include/stdio.h" 3 4
extern FILE *popen (const char *__command, const char *__modes) ;





extern int pclose (FILE *__stream);





extern char *ctermid (char *__s) __attribute__ ((__nothrow__ , __leaf__));
# 913 "/usr/include/stdio.h" 3 4
extern void flockfile (FILE *__stream) __attribute__ ((__nothrow__ , __leaf__));



extern int ftrylockfile (FILE *__stream) __attribute__ ((__nothrow__ , __leaf__)) ;


extern void funlockfile (FILE *__stream) __attribute__ ((__nothrow__ , __leaf__));
# 943 "/usr/include/stdio.h" 3 4

# 139 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h" 2
# 166 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"

# 166 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
typedef int PetscErrorCode;
# 182 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
typedef int PetscClassId;
# 199 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
typedef int PetscMPIInt;
# 208 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
typedef enum { ENUM_DUMMY } PetscEnum;
extern MPI_Datatype MPIU_ENUM ;
# 225 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
typedef int64_t Petsc64bitInt;
# 241 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
typedef int PetscInt;
# 281 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
typedef int PetscBLASInt;
# 292 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
typedef enum { PETSC_PRECISION_SINGLE=4,PETSC_PRECISION_DOUBLE=8 } PetscPrecision;
extern const char *PetscPrecisions[];
# 313 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
extern FILE* PETSC_STDOUT;





extern FILE* PETSC_STDERR;
# 386 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
typedef enum { PETSC_FALSE,PETSC_TRUE } PetscBool;
extern const char *const PetscBools[];
extern MPI_Datatype MPIU_BOOL ;




# 1 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmath.h" 1
# 13 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmath.h"
# 1 "/usr/include/math.h" 1 3 4
# 29 "/usr/include/math.h" 3 4




# 1 "/usr/include/bits/huge_val.h" 1 3 4
# 34 "/usr/include/math.h" 2 3 4

# 1 "/usr/include/bits/huge_valf.h" 1 3 4
# 36 "/usr/include/math.h" 2 3 4
# 1 "/usr/include/bits/huge_vall.h" 1 3 4
# 37 "/usr/include/math.h" 2 3 4


# 1 "/usr/include/bits/inf.h" 1 3 4
# 40 "/usr/include/math.h" 2 3 4


# 1 "/usr/include/bits/nan.h" 1 3 4
# 43 "/usr/include/math.h" 2 3 4



# 1 "/usr/include/bits/mathdef.h" 1 3 4
# 28 "/usr/include/bits/mathdef.h" 3 4

# 28 "/usr/include/bits/mathdef.h" 3 4
typedef float float_t;
typedef double double_t;
# 47 "/usr/include/math.h" 2 3 4
# 70 "/usr/include/math.h" 3 4
# 1 "/usr/include/bits/mathcalls.h" 1 3 4
# 52 "/usr/include/bits/mathcalls.h" 3 4


extern double acos (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __acos (double __x) __attribute__ ((__nothrow__ , __leaf__));

extern double asin (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __asin (double __x) __attribute__ ((__nothrow__ , __leaf__));

extern double atan (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __atan (double __x) __attribute__ ((__nothrow__ , __leaf__));

extern double atan2 (double __y, double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __atan2 (double __y, double __x) __attribute__ ((__nothrow__ , __leaf__));


extern double cos (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __cos (double __x) __attribute__ ((__nothrow__ , __leaf__));

extern double sin (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __sin (double __x) __attribute__ ((__nothrow__ , __leaf__));

extern double tan (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __tan (double __x) __attribute__ ((__nothrow__ , __leaf__));




extern double cosh (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __cosh (double __x) __attribute__ ((__nothrow__ , __leaf__));

extern double sinh (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __sinh (double __x) __attribute__ ((__nothrow__ , __leaf__));

extern double tanh (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __tanh (double __x) __attribute__ ((__nothrow__ , __leaf__));

# 86 "/usr/include/bits/mathcalls.h" 3 4


extern double acosh (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __acosh (double __x) __attribute__ ((__nothrow__ , __leaf__));

extern double asinh (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __asinh (double __x) __attribute__ ((__nothrow__ , __leaf__));

extern double atanh (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __atanh (double __x) __attribute__ ((__nothrow__ , __leaf__));







extern double exp (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __exp (double __x) __attribute__ ((__nothrow__ , __leaf__));


extern double frexp (double __x, int *__exponent) __attribute__ ((__nothrow__ , __leaf__)); extern double __frexp (double __x, int *__exponent) __attribute__ ((__nothrow__ , __leaf__));


extern double ldexp (double __x, int __exponent) __attribute__ ((__nothrow__ , __leaf__)); extern double __ldexp (double __x, int __exponent) __attribute__ ((__nothrow__ , __leaf__));


extern double log (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __log (double __x) __attribute__ ((__nothrow__ , __leaf__));


extern double log10 (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __log10 (double __x) __attribute__ ((__nothrow__ , __leaf__));


extern double modf (double __x, double *__iptr) __attribute__ ((__nothrow__ , __leaf__)); extern double __modf (double __x, double *__iptr) __attribute__ ((__nothrow__ , __leaf__))
     __attribute__ ((__nonnull__ (2)));

# 127 "/usr/include/bits/mathcalls.h" 3 4


extern double expm1 (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __expm1 (double __x) __attribute__ ((__nothrow__ , __leaf__));


extern double log1p (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __log1p (double __x) __attribute__ ((__nothrow__ , __leaf__));


extern double logb (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __logb (double __x) __attribute__ ((__nothrow__ , __leaf__));






extern double exp2 (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __exp2 (double __x) __attribute__ ((__nothrow__ , __leaf__));


extern double log2 (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __log2 (double __x) __attribute__ ((__nothrow__ , __leaf__));








extern double pow (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)); extern double __pow (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__));


extern double sqrt (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __sqrt (double __x) __attribute__ ((__nothrow__ , __leaf__));





extern double hypot (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)); extern double __hypot (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__));






extern double cbrt (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __cbrt (double __x) __attribute__ ((__nothrow__ , __leaf__));








extern double ceil (double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern double __ceil (double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));


extern double fabs (double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern double __fabs (double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));


extern double floor (double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern double __floor (double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));


extern double fmod (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)); extern double __fmod (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__));




extern int __isinf (double __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));


extern int __finite (double __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));





extern int isinf (double __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));


extern int finite (double __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));


extern double drem (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)); extern double __drem (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__));



extern double significand (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __significand (double __x) __attribute__ ((__nothrow__ , __leaf__));





extern double copysign (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern double __copysign (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));






extern double nan (const char *__tagb) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern double __nan (const char *__tagb) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));





extern int __isnan (double __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));



extern int isnan (double __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));


extern double j0 (double) __attribute__ ((__nothrow__ , __leaf__)); extern double __j0 (double) __attribute__ ((__nothrow__ , __leaf__));
extern double j1 (double) __attribute__ ((__nothrow__ , __leaf__)); extern double __j1 (double) __attribute__ ((__nothrow__ , __leaf__));
extern double jn (int, double) __attribute__ ((__nothrow__ , __leaf__)); extern double __jn (int, double) __attribute__ ((__nothrow__ , __leaf__));
extern double y0 (double) __attribute__ ((__nothrow__ , __leaf__)); extern double __y0 (double) __attribute__ ((__nothrow__ , __leaf__));
extern double y1 (double) __attribute__ ((__nothrow__ , __leaf__)); extern double __y1 (double) __attribute__ ((__nothrow__ , __leaf__));
extern double yn (int, double) __attribute__ ((__nothrow__ , __leaf__)); extern double __yn (int, double) __attribute__ ((__nothrow__ , __leaf__));






extern double erf (double) __attribute__ ((__nothrow__ , __leaf__)); extern double __erf (double) __attribute__ ((__nothrow__ , __leaf__));
extern double erfc (double) __attribute__ ((__nothrow__ , __leaf__)); extern double __erfc (double) __attribute__ ((__nothrow__ , __leaf__));
extern double lgamma (double) __attribute__ ((__nothrow__ , __leaf__)); extern double __lgamma (double) __attribute__ ((__nothrow__ , __leaf__));






extern double tgamma (double) __attribute__ ((__nothrow__ , __leaf__)); extern double __tgamma (double) __attribute__ ((__nothrow__ , __leaf__));





extern double gamma (double) __attribute__ ((__nothrow__ , __leaf__)); extern double __gamma (double) __attribute__ ((__nothrow__ , __leaf__));






extern double lgamma_r (double, int *__signgamp) __attribute__ ((__nothrow__ , __leaf__)); extern double __lgamma_r (double, int *__signgamp) __attribute__ ((__nothrow__ , __leaf__));







extern double rint (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __rint (double __x) __attribute__ ((__nothrow__ , __leaf__));


extern double nextafter (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern double __nextafter (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));

extern double nexttoward (double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern double __nexttoward (double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));



extern double remainder (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)); extern double __remainder (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__));



extern double scalbn (double __x, int __n) __attribute__ ((__nothrow__ , __leaf__)); extern double __scalbn (double __x, int __n) __attribute__ ((__nothrow__ , __leaf__));



extern int ilogb (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern int __ilogb (double __x) __attribute__ ((__nothrow__ , __leaf__));




extern double scalbln (double __x, long int __n) __attribute__ ((__nothrow__ , __leaf__)); extern double __scalbln (double __x, long int __n) __attribute__ ((__nothrow__ , __leaf__));



extern double nearbyint (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __nearbyint (double __x) __attribute__ ((__nothrow__ , __leaf__));



extern double round (double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern double __round (double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));



extern double trunc (double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern double __trunc (double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));




extern double remquo (double __x, double __y, int *__quo) __attribute__ ((__nothrow__ , __leaf__)); extern double __remquo (double __x, double __y, int *__quo) __attribute__ ((__nothrow__ , __leaf__));






extern long int lrint (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long int __lrint (double __x) __attribute__ ((__nothrow__ , __leaf__));
extern long long int llrint (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long long int __llrint (double __x) __attribute__ ((__nothrow__ , __leaf__));



extern long int lround (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long int __lround (double __x) __attribute__ ((__nothrow__ , __leaf__));
extern long long int llround (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long long int __llround (double __x) __attribute__ ((__nothrow__ , __leaf__));



extern double fdim (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)); extern double __fdim (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__));


extern double fmax (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern double __fmax (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));


extern double fmin (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern double __fmin (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));



extern int __fpclassify (double __value) __attribute__ ((__nothrow__ , __leaf__))
     __attribute__ ((__const__));


extern int __signbit (double __value) __attribute__ ((__nothrow__ , __leaf__))
     __attribute__ ((__const__));



extern double fma (double __x, double __y, double __z) __attribute__ ((__nothrow__ , __leaf__)); extern double __fma (double __x, double __y, double __z) __attribute__ ((__nothrow__ , __leaf__));








extern double scalb (double __x, double __n) __attribute__ ((__nothrow__ , __leaf__)); extern double __scalb (double __x, double __n) __attribute__ ((__nothrow__ , __leaf__));
# 71 "/usr/include/math.h" 2 3 4
# 89 "/usr/include/math.h" 3 4
# 1 "/usr/include/bits/mathcalls.h" 1 3 4
# 52 "/usr/include/bits/mathcalls.h" 3 4


extern float acosf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __acosf (float __x) __attribute__ ((__nothrow__ , __leaf__));

extern float asinf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __asinf (float __x) __attribute__ ((__nothrow__ , __leaf__));

extern float atanf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __atanf (float __x) __attribute__ ((__nothrow__ , __leaf__));

extern float atan2f (float __y, float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __atan2f (float __y, float __x) __attribute__ ((__nothrow__ , __leaf__));


extern float cosf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __cosf (float __x) __attribute__ ((__nothrow__ , __leaf__));

extern float sinf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __sinf (float __x) __attribute__ ((__nothrow__ , __leaf__));

extern float tanf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __tanf (float __x) __attribute__ ((__nothrow__ , __leaf__));




extern float coshf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __coshf (float __x) __attribute__ ((__nothrow__ , __leaf__));

extern float sinhf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __sinhf (float __x) __attribute__ ((__nothrow__ , __leaf__));

extern float tanhf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __tanhf (float __x) __attribute__ ((__nothrow__ , __leaf__));

# 86 "/usr/include/bits/mathcalls.h" 3 4


extern float acoshf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __acoshf (float __x) __attribute__ ((__nothrow__ , __leaf__));

extern float asinhf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __asinhf (float __x) __attribute__ ((__nothrow__ , __leaf__));

extern float atanhf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __atanhf (float __x) __attribute__ ((__nothrow__ , __leaf__));







extern float expf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __expf (float __x) __attribute__ ((__nothrow__ , __leaf__));


extern float frexpf (float __x, int *__exponent) __attribute__ ((__nothrow__ , __leaf__)); extern float __frexpf (float __x, int *__exponent) __attribute__ ((__nothrow__ , __leaf__));


extern float ldexpf (float __x, int __exponent) __attribute__ ((__nothrow__ , __leaf__)); extern float __ldexpf (float __x, int __exponent) __attribute__ ((__nothrow__ , __leaf__));


extern float logf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __logf (float __x) __attribute__ ((__nothrow__ , __leaf__));


extern float log10f (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __log10f (float __x) __attribute__ ((__nothrow__ , __leaf__));


extern float modff (float __x, float *__iptr) __attribute__ ((__nothrow__ , __leaf__)); extern float __modff (float __x, float *__iptr) __attribute__ ((__nothrow__ , __leaf__))
     __attribute__ ((__nonnull__ (2)));

# 127 "/usr/include/bits/mathcalls.h" 3 4


extern float expm1f (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __expm1f (float __x) __attribute__ ((__nothrow__ , __leaf__));


extern float log1pf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __log1pf (float __x) __attribute__ ((__nothrow__ , __leaf__));


extern float logbf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __logbf (float __x) __attribute__ ((__nothrow__ , __leaf__));






extern float exp2f (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __exp2f (float __x) __attribute__ ((__nothrow__ , __leaf__));


extern float log2f (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __log2f (float __x) __attribute__ ((__nothrow__ , __leaf__));








extern float powf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)); extern float __powf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__));


extern float sqrtf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __sqrtf (float __x) __attribute__ ((__nothrow__ , __leaf__));





extern float hypotf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)); extern float __hypotf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__));






extern float cbrtf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __cbrtf (float __x) __attribute__ ((__nothrow__ , __leaf__));








extern float ceilf (float __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern float __ceilf (float __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));


extern float fabsf (float __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern float __fabsf (float __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));


extern float floorf (float __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern float __floorf (float __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));


extern float fmodf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)); extern float __fmodf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__));




extern int __isinff (float __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));


extern int __finitef (float __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));





extern int isinff (float __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));


extern int finitef (float __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));


extern float dremf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)); extern float __dremf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__));



extern float significandf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __significandf (float __x) __attribute__ ((__nothrow__ , __leaf__));





extern float copysignf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern float __copysignf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));






extern float nanf (const char *__tagb) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern float __nanf (const char *__tagb) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));





extern int __isnanf (float __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));



extern int isnanf (float __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));


extern float j0f (float) __attribute__ ((__nothrow__ , __leaf__)); extern float __j0f (float) __attribute__ ((__nothrow__ , __leaf__));
extern float j1f (float) __attribute__ ((__nothrow__ , __leaf__)); extern float __j1f (float) __attribute__ ((__nothrow__ , __leaf__));
extern float jnf (int, float) __attribute__ ((__nothrow__ , __leaf__)); extern float __jnf (int, float) __attribute__ ((__nothrow__ , __leaf__));
extern float y0f (float) __attribute__ ((__nothrow__ , __leaf__)); extern float __y0f (float) __attribute__ ((__nothrow__ , __leaf__));
extern float y1f (float) __attribute__ ((__nothrow__ , __leaf__)); extern float __y1f (float) __attribute__ ((__nothrow__ , __leaf__));
extern float ynf (int, float) __attribute__ ((__nothrow__ , __leaf__)); extern float __ynf (int, float) __attribute__ ((__nothrow__ , __leaf__));






extern float erff (float) __attribute__ ((__nothrow__ , __leaf__)); extern float __erff (float) __attribute__ ((__nothrow__ , __leaf__));
extern float erfcf (float) __attribute__ ((__nothrow__ , __leaf__)); extern float __erfcf (float) __attribute__ ((__nothrow__ , __leaf__));
extern float lgammaf (float) __attribute__ ((__nothrow__ , __leaf__)); extern float __lgammaf (float) __attribute__ ((__nothrow__ , __leaf__));






extern float tgammaf (float) __attribute__ ((__nothrow__ , __leaf__)); extern float __tgammaf (float) __attribute__ ((__nothrow__ , __leaf__));





extern float gammaf (float) __attribute__ ((__nothrow__ , __leaf__)); extern float __gammaf (float) __attribute__ ((__nothrow__ , __leaf__));






extern float lgammaf_r (float, int *__signgamp) __attribute__ ((__nothrow__ , __leaf__)); extern float __lgammaf_r (float, int *__signgamp) __attribute__ ((__nothrow__ , __leaf__));







extern float rintf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __rintf (float __x) __attribute__ ((__nothrow__ , __leaf__));


extern float nextafterf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern float __nextafterf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));

extern float nexttowardf (float __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern float __nexttowardf (float __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));



extern float remainderf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)); extern float __remainderf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__));



extern float scalbnf (float __x, int __n) __attribute__ ((__nothrow__ , __leaf__)); extern float __scalbnf (float __x, int __n) __attribute__ ((__nothrow__ , __leaf__));



extern int ilogbf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern int __ilogbf (float __x) __attribute__ ((__nothrow__ , __leaf__));




extern float scalblnf (float __x, long int __n) __attribute__ ((__nothrow__ , __leaf__)); extern float __scalblnf (float __x, long int __n) __attribute__ ((__nothrow__ , __leaf__));



extern float nearbyintf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __nearbyintf (float __x) __attribute__ ((__nothrow__ , __leaf__));



extern float roundf (float __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern float __roundf (float __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));



extern float truncf (float __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern float __truncf (float __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));




extern float remquof (float __x, float __y, int *__quo) __attribute__ ((__nothrow__ , __leaf__)); extern float __remquof (float __x, float __y, int *__quo) __attribute__ ((__nothrow__ , __leaf__));






extern long int lrintf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern long int __lrintf (float __x) __attribute__ ((__nothrow__ , __leaf__));
extern long long int llrintf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern long long int __llrintf (float __x) __attribute__ ((__nothrow__ , __leaf__));



extern long int lroundf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern long int __lroundf (float __x) __attribute__ ((__nothrow__ , __leaf__));
extern long long int llroundf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern long long int __llroundf (float __x) __attribute__ ((__nothrow__ , __leaf__));



extern float fdimf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)); extern float __fdimf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__));


extern float fmaxf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern float __fmaxf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));


extern float fminf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern float __fminf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));



extern int __fpclassifyf (float __value) __attribute__ ((__nothrow__ , __leaf__))
     __attribute__ ((__const__));


extern int __signbitf (float __value) __attribute__ ((__nothrow__ , __leaf__))
     __attribute__ ((__const__));



extern float fmaf (float __x, float __y, float __z) __attribute__ ((__nothrow__ , __leaf__)); extern float __fmaf (float __x, float __y, float __z) __attribute__ ((__nothrow__ , __leaf__));








extern float scalbf (float __x, float __n) __attribute__ ((__nothrow__ , __leaf__)); extern float __scalbf (float __x, float __n) __attribute__ ((__nothrow__ , __leaf__));
# 90 "/usr/include/math.h" 2 3 4
# 133 "/usr/include/math.h" 3 4
# 1 "/usr/include/bits/mathcalls.h" 1 3 4
# 52 "/usr/include/bits/mathcalls.h" 3 4


extern long double acosl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __acosl (long double __x) __attribute__ ((__nothrow__ , __leaf__));

extern long double asinl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __asinl (long double __x) __attribute__ ((__nothrow__ , __leaf__));

extern long double atanl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __atanl (long double __x) __attribute__ ((__nothrow__ , __leaf__));

extern long double atan2l (long double __y, long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __atan2l (long double __y, long double __x) __attribute__ ((__nothrow__ , __leaf__));


extern long double cosl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __cosl (long double __x) __attribute__ ((__nothrow__ , __leaf__));

extern long double sinl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __sinl (long double __x) __attribute__ ((__nothrow__ , __leaf__));

extern long double tanl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __tanl (long double __x) __attribute__ ((__nothrow__ , __leaf__));




extern long double coshl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __coshl (long double __x) __attribute__ ((__nothrow__ , __leaf__));

extern long double sinhl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __sinhl (long double __x) __attribute__ ((__nothrow__ , __leaf__));

extern long double tanhl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __tanhl (long double __x) __attribute__ ((__nothrow__ , __leaf__));

# 86 "/usr/include/bits/mathcalls.h" 3 4


extern long double acoshl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __acoshl (long double __x) __attribute__ ((__nothrow__ , __leaf__));

extern long double asinhl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __asinhl (long double __x) __attribute__ ((__nothrow__ , __leaf__));

extern long double atanhl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __atanhl (long double __x) __attribute__ ((__nothrow__ , __leaf__));







extern long double expl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __expl (long double __x) __attribute__ ((__nothrow__ , __leaf__));


extern long double frexpl (long double __x, int *__exponent) __attribute__ ((__nothrow__ , __leaf__)); extern long double __frexpl (long double __x, int *__exponent) __attribute__ ((__nothrow__ , __leaf__));


extern long double ldexpl (long double __x, int __exponent) __attribute__ ((__nothrow__ , __leaf__)); extern long double __ldexpl (long double __x, int __exponent) __attribute__ ((__nothrow__ , __leaf__));


extern long double logl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __logl (long double __x) __attribute__ ((__nothrow__ , __leaf__));


extern long double log10l (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __log10l (long double __x) __attribute__ ((__nothrow__ , __leaf__));


extern long double modfl (long double __x, long double *__iptr) __attribute__ ((__nothrow__ , __leaf__)); extern long double __modfl (long double __x, long double *__iptr) __attribute__ ((__nothrow__ , __leaf__))
     __attribute__ ((__nonnull__ (2)));

# 127 "/usr/include/bits/mathcalls.h" 3 4


extern long double expm1l (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __expm1l (long double __x) __attribute__ ((__nothrow__ , __leaf__));


extern long double log1pl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __log1pl (long double __x) __attribute__ ((__nothrow__ , __leaf__));


extern long double logbl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __logbl (long double __x) __attribute__ ((__nothrow__ , __leaf__));






extern long double exp2l (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __exp2l (long double __x) __attribute__ ((__nothrow__ , __leaf__));


extern long double log2l (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __log2l (long double __x) __attribute__ ((__nothrow__ , __leaf__));








extern long double powl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)); extern long double __powl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__));


extern long double sqrtl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __sqrtl (long double __x) __attribute__ ((__nothrow__ , __leaf__));





extern long double hypotl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)); extern long double __hypotl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__));






extern long double cbrtl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __cbrtl (long double __x) __attribute__ ((__nothrow__ , __leaf__));








extern long double ceill (long double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern long double __ceill (long double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));


extern long double fabsl (long double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern long double __fabsl (long double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));


extern long double floorl (long double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern long double __floorl (long double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));


extern long double fmodl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)); extern long double __fmodl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__));




extern int __isinfl (long double __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));


extern int __finitel (long double __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));





extern int isinfl (long double __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));


extern int finitel (long double __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));


extern long double dreml (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)); extern long double __dreml (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__));



extern long double significandl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __significandl (long double __x) __attribute__ ((__nothrow__ , __leaf__));





extern long double copysignl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern long double __copysignl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));






extern long double nanl (const char *__tagb) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern long double __nanl (const char *__tagb) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));





extern int __isnanl (long double __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));



extern int isnanl (long double __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));


extern long double j0l (long double) __attribute__ ((__nothrow__ , __leaf__)); extern long double __j0l (long double) __attribute__ ((__nothrow__ , __leaf__));
extern long double j1l (long double) __attribute__ ((__nothrow__ , __leaf__)); extern long double __j1l (long double) __attribute__ ((__nothrow__ , __leaf__));
extern long double jnl (int, long double) __attribute__ ((__nothrow__ , __leaf__)); extern long double __jnl (int, long double) __attribute__ ((__nothrow__ , __leaf__));
extern long double y0l (long double) __attribute__ ((__nothrow__ , __leaf__)); extern long double __y0l (long double) __attribute__ ((__nothrow__ , __leaf__));
extern long double y1l (long double) __attribute__ ((__nothrow__ , __leaf__)); extern long double __y1l (long double) __attribute__ ((__nothrow__ , __leaf__));
extern long double ynl (int, long double) __attribute__ ((__nothrow__ , __leaf__)); extern long double __ynl (int, long double) __attribute__ ((__nothrow__ , __leaf__));






extern long double erfl (long double) __attribute__ ((__nothrow__ , __leaf__)); extern long double __erfl (long double) __attribute__ ((__nothrow__ , __leaf__));
extern long double erfcl (long double) __attribute__ ((__nothrow__ , __leaf__)); extern long double __erfcl (long double) __attribute__ ((__nothrow__ , __leaf__));
extern long double lgammal (long double) __attribute__ ((__nothrow__ , __leaf__)); extern long double __lgammal (long double) __attribute__ ((__nothrow__ , __leaf__));






extern long double tgammal (long double) __attribute__ ((__nothrow__ , __leaf__)); extern long double __tgammal (long double) __attribute__ ((__nothrow__ , __leaf__));





extern long double gammal (long double) __attribute__ ((__nothrow__ , __leaf__)); extern long double __gammal (long double) __attribute__ ((__nothrow__ , __leaf__));






extern long double lgammal_r (long double, int *__signgamp) __attribute__ ((__nothrow__ , __leaf__)); extern long double __lgammal_r (long double, int *__signgamp) __attribute__ ((__nothrow__ , __leaf__));







extern long double rintl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __rintl (long double __x) __attribute__ ((__nothrow__ , __leaf__));


extern long double nextafterl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern long double __nextafterl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));

extern long double nexttowardl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern long double __nexttowardl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));



extern long double remainderl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)); extern long double __remainderl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__));



extern long double scalbnl (long double __x, int __n) __attribute__ ((__nothrow__ , __leaf__)); extern long double __scalbnl (long double __x, int __n) __attribute__ ((__nothrow__ , __leaf__));



extern int ilogbl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern int __ilogbl (long double __x) __attribute__ ((__nothrow__ , __leaf__));




extern long double scalblnl (long double __x, long int __n) __attribute__ ((__nothrow__ , __leaf__)); extern long double __scalblnl (long double __x, long int __n) __attribute__ ((__nothrow__ , __leaf__));



extern long double nearbyintl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __nearbyintl (long double __x) __attribute__ ((__nothrow__ , __leaf__));



extern long double roundl (long double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern long double __roundl (long double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));



extern long double truncl (long double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern long double __truncl (long double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));




extern long double remquol (long double __x, long double __y, int *__quo) __attribute__ ((__nothrow__ , __leaf__)); extern long double __remquol (long double __x, long double __y, int *__quo) __attribute__ ((__nothrow__ , __leaf__));






extern long int lrintl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long int __lrintl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
extern long long int llrintl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long long int __llrintl (long double __x) __attribute__ ((__nothrow__ , __leaf__));



extern long int lroundl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long int __lroundl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
extern long long int llroundl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long long int __llroundl (long double __x) __attribute__ ((__nothrow__ , __leaf__));



extern long double fdiml (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)); extern long double __fdiml (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__));


extern long double fmaxl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern long double __fmaxl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));


extern long double fminl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern long double __fminl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));



extern int __fpclassifyl (long double __value) __attribute__ ((__nothrow__ , __leaf__))
     __attribute__ ((__const__));


extern int __signbitl (long double __value) __attribute__ ((__nothrow__ , __leaf__))
     __attribute__ ((__const__));



extern long double fmal (long double __x, long double __y, long double __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double __fmal (long double __x, long double __y, long double __z) __attribute__ ((__nothrow__ , __leaf__));








extern long double scalbl (long double __x, long double __n) __attribute__ ((__nothrow__ , __leaf__)); extern long double __scalbl (long double __x, long double __n) __attribute__ ((__nothrow__ , __leaf__));
# 134 "/usr/include/math.h" 2 3 4
# 149 "/usr/include/math.h" 3 4
extern int signgam;
# 190 "/usr/include/math.h" 3 4
enum
  {
    FP_NAN =

      0,
    FP_INFINITE =

      1,
    FP_ZERO =

      2,
    FP_SUBNORMAL =

      3,
    FP_NORMAL =

      4
  };
# 288 "/usr/include/math.h" 3 4
typedef enum
{
  _IEEE_ = -1,
  _SVID_,
  _XOPEN_,
  _POSIX_,
  _ISOC_
} _LIB_VERSION_TYPE;




extern _LIB_VERSION_TYPE _LIB_VERSION;
# 313 "/usr/include/math.h" 3 4
struct exception

  {
    int type;
    char *name;
    double arg1;
    double arg2;
    double retval;
  };




extern int matherr (struct exception *__exc);
# 475 "/usr/include/math.h" 3 4

# 14 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmath.h" 2
# 52 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmath.h"

# 52 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmath.h"
typedef double PetscReal;
# 155 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmath.h"
# 1 "/usr/include/complex.h" 1 3 4
# 28 "/usr/include/complex.h" 3 4
# 1 "/usr/include/bits/mathdef.h" 1 3 4
# 29 "/usr/include/complex.h" 2 3 4


# 75 "/usr/include/complex.h" 3 4
# 1 "/usr/include/bits/cmathcalls.h" 1 3 4
# 53 "/usr/include/bits/cmathcalls.h" 3 4

# 53 "/usr/include/bits/cmathcalls.h" 3 4
extern double _Complex cacos (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double _Complex __cacos (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

extern double _Complex casin (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double _Complex __casin (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

extern double _Complex catan (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double _Complex __catan (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));


extern double _Complex ccos (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double _Complex __ccos (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

extern double _Complex csin (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double _Complex __csin (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

extern double _Complex ctan (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double _Complex __ctan (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));





extern double _Complex cacosh (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double _Complex __cacosh (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

extern double _Complex casinh (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double _Complex __casinh (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

extern double _Complex catanh (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double _Complex __catanh (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));


extern double _Complex ccosh (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double _Complex __ccosh (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

extern double _Complex csinh (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double _Complex __csinh (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

extern double _Complex ctanh (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double _Complex __ctanh (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));





extern double _Complex cexp (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double _Complex __cexp (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));


extern double _Complex clog (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double _Complex __clog (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));
# 101 "/usr/include/bits/cmathcalls.h" 3 4
extern double _Complex cpow (double _Complex __x, double _Complex __y) __attribute__ ((__nothrow__ , __leaf__)); extern double _Complex __cpow (double _Complex __x, double _Complex __y) __attribute__ ((__nothrow__ , __leaf__));


extern double _Complex csqrt (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double _Complex __csqrt (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));





extern double cabs (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double __cabs (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));


extern double carg (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double __carg (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));


extern double _Complex conj (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double _Complex __conj (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));


extern double _Complex cproj (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double _Complex __cproj (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));





extern double cimag (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double __cimag (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));


extern double creal (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double __creal (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));
# 76 "/usr/include/complex.h" 2 3 4
# 85 "/usr/include/complex.h" 3 4
# 1 "/usr/include/bits/cmathcalls.h" 1 3 4
# 53 "/usr/include/bits/cmathcalls.h" 3 4
extern float _Complex cacosf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float _Complex __cacosf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

extern float _Complex casinf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float _Complex __casinf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

extern float _Complex catanf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float _Complex __catanf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));


extern float _Complex ccosf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float _Complex __ccosf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

extern float _Complex csinf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float _Complex __csinf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

extern float _Complex ctanf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float _Complex __ctanf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));





extern float _Complex cacoshf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float _Complex __cacoshf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

extern float _Complex casinhf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float _Complex __casinhf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

extern float _Complex catanhf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float _Complex __catanhf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));


extern float _Complex ccoshf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float _Complex __ccoshf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

extern float _Complex csinhf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float _Complex __csinhf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

extern float _Complex ctanhf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float _Complex __ctanhf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));





extern float _Complex cexpf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float _Complex __cexpf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));


extern float _Complex clogf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float _Complex __clogf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));
# 101 "/usr/include/bits/cmathcalls.h" 3 4
extern float _Complex cpowf (float _Complex __x, float _Complex __y) __attribute__ ((__nothrow__ , __leaf__)); extern float _Complex __cpowf (float _Complex __x, float _Complex __y) __attribute__ ((__nothrow__ , __leaf__));


extern float _Complex csqrtf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float _Complex __csqrtf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));





extern float cabsf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float __cabsf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));


extern float cargf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float __cargf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));


extern float _Complex conjf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float _Complex __conjf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));


extern float _Complex cprojf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float _Complex __cprojf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));





extern float cimagf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float __cimagf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));


extern float crealf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float __crealf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));
# 86 "/usr/include/complex.h" 2 3 4
# 104 "/usr/include/complex.h" 3 4
# 1 "/usr/include/bits/cmathcalls.h" 1 3 4
# 53 "/usr/include/bits/cmathcalls.h" 3 4
extern long double _Complex cacosl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double _Complex __cacosl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

extern long double _Complex casinl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double _Complex __casinl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

extern long double _Complex catanl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double _Complex __catanl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));


extern long double _Complex ccosl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double _Complex __ccosl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

extern long double _Complex csinl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double _Complex __csinl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

extern long double _Complex ctanl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double _Complex __ctanl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));





extern long double _Complex cacoshl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double _Complex __cacoshl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

extern long double _Complex casinhl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double _Complex __casinhl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

extern long double _Complex catanhl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double _Complex __catanhl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));


extern long double _Complex ccoshl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double _Complex __ccoshl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

extern long double _Complex csinhl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double _Complex __csinhl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

extern long double _Complex ctanhl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double _Complex __ctanhl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));





extern long double _Complex cexpl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double _Complex __cexpl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));


extern long double _Complex clogl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double _Complex __clogl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));
# 101 "/usr/include/bits/cmathcalls.h" 3 4
extern long double _Complex cpowl (long double _Complex __x, long double _Complex __y) __attribute__ ((__nothrow__ , __leaf__)); extern long double _Complex __cpowl (long double _Complex __x, long double _Complex __y) __attribute__ ((__nothrow__ , __leaf__));


extern long double _Complex csqrtl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double _Complex __csqrtl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));





extern long double cabsl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double __cabsl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));


extern long double cargl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double __cargl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));


extern long double _Complex conjl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double _Complex __conjl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));


extern long double _Complex cprojl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double _Complex __cprojl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));





extern long double cimagl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double __cimagl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));


extern long double creall (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double __creall (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));
# 105 "/usr/include/complex.h" 2 3 4








# 156 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmath.h" 2
# 178 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmath.h"

# 178 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmath.h"
typedef double _Complex PetscComplex;
# 278 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmath.h"
typedef PetscReal PetscScalar;




static inline PetscReal PetscAbsScalar(PetscScalar a) {return a < 0.0 ? -a : a;}
# 324 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmath.h"
typedef enum { PETSC_SCALAR_DOUBLE,PETSC_SCALAR_SINGLE, PETSC_SCALAR_LONG_DOUBLE } PetscScalarPrecision;



extern PetscComplex PETSC_i;
# 501 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmath.h"
extern PetscErrorCode PetscIsInfOrNanReal(PetscReal);
extern PetscBool PetscIsNormalReal(PetscReal);
static inline PetscErrorCode PetscIsInfOrNanScalar(PetscScalar v) {return PetscIsInfOrNanReal(PetscAbsScalar(v));}
static inline PetscErrorCode PetscIsNormalScalar(PetscScalar v) {return PetscIsNormalReal(PetscAbsScalar(v));}
# 515 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmath.h"
typedef PetscScalar MatScalar;
typedef PetscReal MatReal;

struct petsc_mpiu_2scalar {PetscScalar a,b;};
extern MPI_Datatype MPIU_2SCALAR ;







static inline PetscInt PetscPowInt(PetscInt base,PetscInt power)
{
  PetscInt result = 1;
  while (power) {
    if (power & 1) result *= base;
    power >>= 1;
    base *= base;
  }
  return result;
}

static inline PetscReal PetscPowRealInt(PetscReal base,PetscInt power)
{
  PetscReal result = 1;
  if (power < 0) {
    power = -power;
    if (base != 0.0) base = 1./base;
  }
  while (power) {
    if (power & 1) result *= base;
    power >>= 1;
    base *= base;
  }
  return result;
}

static inline PetscScalar PetscPowScalarInt(PetscScalar base,PetscInt power)
{
  PetscScalar result = 1;
  if (power < 0) {
    power = -power;
    if (base != 0.0) base = 1./base;
  }
  while (power) {
    if (power & 1) result *= base;
    power >>= 1;
    base *= base;
  }
  return result;
}

static inline PetscScalar PetscPowScalarReal(PetscScalar base,PetscReal power)
{
  PetscScalar cpower = power;
  return pow(base,cpower);
}
# 394 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h" 2
# 407 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
typedef enum { PETSC_COPY_VALUES, PETSC_OWN_POINTER, PETSC_USE_POINTER} PetscCopyMode;
extern const char *const PetscCopyModes[];
# 518 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
extern MPI_Comm PETSC_COMM_WORLD;
# 532 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
extern PetscBool PetscBeganMPI;
extern PetscBool PetscInitializeCalled;
extern PetscBool PetscFinalizeCalled;
extern PetscBool PetscCUSPSynchronize;
extern PetscBool PetscViennaCLSynchronize;

extern PetscErrorCode PetscSetHelpVersionFunctions(PetscErrorCode (*)(MPI_Comm),PetscErrorCode (*)(MPI_Comm));
extern PetscErrorCode PetscCommDuplicate(MPI_Comm,MPI_Comm*,int*);
extern PetscErrorCode PetscCommDestroy(MPI_Comm*);
# 1350 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
extern PetscErrorCode (*PetscTrMalloc)(size_t,int,const char[],const char[],void**);
extern PetscErrorCode (*PetscTrFree)(void*,int,const char[],const char[]);
extern PetscErrorCode PetscMallocSet(PetscErrorCode (*)(size_t,int,const char[],const char[],void**),PetscErrorCode (*)(void*,int,const char[],const char[]));
extern PetscErrorCode PetscMallocClear(void);






typedef double PetscLogDouble;





extern PetscErrorCode PetscMallocDump(FILE *);
extern PetscErrorCode PetscMallocDumpLog(FILE *);
extern PetscErrorCode PetscMallocGetCurrentUsage(PetscLogDouble *);
extern PetscErrorCode PetscMallocGetMaximumUsage(PetscLogDouble *);
extern PetscErrorCode PetscMallocDebug(PetscBool);
extern PetscErrorCode PetscMallocGetDebug(PetscBool*);
extern PetscErrorCode PetscMallocValidate(int,const char[],const char[]);
extern PetscErrorCode PetscMallocSetDumpLog(void);
extern PetscErrorCode PetscMallocSetDumpLogThreshold(PetscLogDouble);
extern PetscErrorCode PetscMallocGetDumpLog(PetscBool*);
# 1388 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
typedef enum {PETSC_INT = 0,PETSC_DOUBLE = 1,PETSC_COMPLEX = 2, PETSC_LONG = 3 ,PETSC_SHORT = 4,PETSC_FLOAT = 5,
              PETSC_CHAR = 6,PETSC_BIT_LOGICAL = 7,PETSC_ENUM = 8,PETSC_BOOL=9, PETSC___FLOAT128 = 10,PETSC_OBJECT = 11, PETSC_FUNCTION = 12, PETSC_STRING = 12} PetscDataType;
extern const char *const PetscDataTypes[];
# 1412 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
extern PetscErrorCode PetscDataTypeToMPIDataType(PetscDataType,MPI_Datatype*);
extern PetscErrorCode PetscMPIDataTypeToPetscDataType(MPI_Datatype,PetscDataType*);
extern PetscErrorCode PetscDataTypeGetSize(PetscDataType,size_t*);
extern PetscErrorCode PetscDataTypeFromString(const char*,PetscDataType*,PetscBool*);






extern PetscErrorCode PetscBitMemcpy(void*,PetscInt,const void*,PetscInt,PetscInt,PetscDataType);
extern PetscErrorCode PetscMemmove(void*,void *,size_t);
extern PetscErrorCode PetscMemcmp(const void*,const void*,size_t,PetscBool *);
extern PetscErrorCode PetscStrlen(const char[],size_t*);
extern PetscErrorCode PetscStrToArray(const char[],char,int*,char ***);
extern PetscErrorCode PetscStrToArrayDestroy(int,char **);
extern PetscErrorCode PetscStrcmp(const char[],const char[],PetscBool *);
extern PetscErrorCode PetscStrgrt(const char[],const char[],PetscBool *);
extern PetscErrorCode PetscStrcasecmp(const char[],const char[],PetscBool *);
extern PetscErrorCode PetscStrncmp(const char[],const char[],size_t,PetscBool *);
extern PetscErrorCode PetscStrcpy(char[],const char[]);
extern PetscErrorCode PetscStrcat(char[],const char[]);
extern PetscErrorCode PetscStrncat(char[],const char[],size_t);
extern PetscErrorCode PetscStrncpy(char[],const char[],size_t);
extern PetscErrorCode PetscStrchr(const char[],char,char *[]);
extern PetscErrorCode PetscStrtolower(char[]);
extern PetscErrorCode PetscStrtoupper(char[]);
extern PetscErrorCode PetscStrrchr(const char[],char,char *[]);
extern PetscErrorCode PetscStrstr(const char[],const char[],char *[]);
extern PetscErrorCode PetscStrrstr(const char[],const char[],char *[]);
extern PetscErrorCode PetscStrendswith(const char[],const char[],PetscBool*);
extern PetscErrorCode PetscStrbeginswith(const char[],const char[],PetscBool*);
extern PetscErrorCode PetscStrendswithwhich(const char[],const char *const*,PetscInt*);
extern PetscErrorCode PetscStrallocpy(const char[],char *[]);
extern PetscErrorCode PetscStrArrayallocpy(const char *const*,char***);
extern PetscErrorCode PetscStrArrayDestroy(char***);
extern PetscErrorCode PetscStrNArrayallocpy(PetscInt,const char *const*,char***);
extern PetscErrorCode PetscStrNArrayDestroy(PetscInt,char***);
extern PetscErrorCode PetscStrreplace(MPI_Comm,const char[],char[],size_t);

extern void PetscStrcmpNoError(const char[],const char[],PetscBool *);
# 1461 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
typedef struct _p_PetscToken* PetscToken;

extern PetscErrorCode PetscTokenCreate(const char[],const char,PetscToken*);
extern PetscErrorCode PetscTokenFind(PetscToken,char *[]);
extern PetscErrorCode PetscTokenDestroy(PetscToken*);

extern PetscErrorCode PetscEListFind(PetscInt,const char *const*,const char*,PetscInt*,PetscBool*);
extern PetscErrorCode PetscEnumFind(const char *const*,const char*,PetscEnum*,PetscBool*);




extern MPI_Op PetscMaxSum_Op;
# 1486 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
extern PetscErrorCode PetscMaxSum(MPI_Comm,const PetscInt[],PetscInt*,PetscInt*);

extern PetscErrorCode MPIULong_Send(void*,PetscInt,MPI_Datatype,PetscMPIInt,PetscMPIInt,MPI_Comm);
extern PetscErrorCode MPIULong_Recv(void*,PetscInt,MPI_Datatype,PetscMPIInt,PetscMPIInt,MPI_Comm);
# 1500 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
typedef struct _p_PetscObject* PetscObject;
# 1511 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
typedef Petsc64bitInt PetscObjectId;
# 1524 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
typedef Petsc64bitInt PetscObjectState;
# 1534 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
typedef struct _n_PetscFunctionList *PetscFunctionList;
# 1553 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
typedef enum {FILE_MODE_READ, FILE_MODE_WRITE, FILE_MODE_APPEND, FILE_MODE_UPDATE, FILE_MODE_APPEND_UPDATE} PetscFileMode;
extern const char *const PetscFileModes[];




# 1 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscerror.h" 1
# 540 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscerror.h"
typedef enum {PETSC_ERROR_INITIAL=0,PETSC_ERROR_REPEAT=1,PETSC_ERROR_IN_CXX = 2} PetscErrorType;




extern PetscErrorCode PetscError(MPI_Comm,int,const char*,const char*,PetscErrorCode,PetscErrorType,const char*,...);

extern PetscErrorCode PetscErrorPrintfInitialize(void);
extern PetscErrorCode PetscErrorMessage(int,const char*[],char **);
extern PetscErrorCode PetscTraceBackErrorHandler(MPI_Comm,int,const char*,const char*,PetscErrorCode,PetscErrorType,const char*,void*);
extern PetscErrorCode PetscIgnoreErrorHandler(MPI_Comm,int,const char*,const char*,PetscErrorCode,PetscErrorType,const char*,void*);
extern PetscErrorCode PetscEmacsClientErrorHandler(MPI_Comm,int,const char*,const char*,PetscErrorCode,PetscErrorType,const char*,void*);
extern PetscErrorCode PetscMPIAbortErrorHandler(MPI_Comm,int,const char*,const char*,PetscErrorCode,PetscErrorType,const char*,void*);
extern PetscErrorCode PetscAbortErrorHandler(MPI_Comm,int,const char*,const char*,PetscErrorCode,PetscErrorType,const char*,void*);
extern PetscErrorCode PetscAttachDebuggerErrorHandler(MPI_Comm,int,const char*,const char*,PetscErrorCode,PetscErrorType,const char*,void*);
extern PetscErrorCode PetscReturnErrorHandler(MPI_Comm,int,const char*,const char*,PetscErrorCode,PetscErrorType,const char*,void*);
extern PetscErrorCode PetscPushErrorHandler(PetscErrorCode (*handler)(MPI_Comm,int,const char*,const char*,PetscErrorCode,PetscErrorType,const char*,void*),void*);
extern PetscErrorCode PetscPopErrorHandler(void);
extern PetscErrorCode PetscSignalHandlerDefault(int,void*);
extern PetscErrorCode PetscPushSignalHandler(PetscErrorCode (*)(int,void *),void*);
extern PetscErrorCode PetscPopSignalHandler(void);
extern PetscErrorCode PetscCheckPointerSetIntensity(PetscInt);
# 601 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscerror.h"
extern PetscErrorCode (*PetscErrorPrintf)(const char[],...);

typedef enum {PETSC_FP_TRAP_OFF=0,PETSC_FP_TRAP_ON=1} PetscFPTrap;
extern PetscErrorCode PetscSetFPTrap(PetscFPTrap);
extern PetscErrorCode PetscFPTrapPush(PetscFPTrap);
extern PetscErrorCode PetscFPTrapPop(void);







typedef struct {
  const char *function[64];
  const char *file[64];
        int line[64];
        PetscBool petscroutine[64];
        int currentsize;
        int hotdepth;
} PetscStack;

extern PetscStack *petscstack;

PetscErrorCode PetscStackCopy(PetscStack*,PetscStack*);
PetscErrorCode PetscStackPrint(PetscStack *,FILE*);
# 846 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscerror.h"
static inline PetscBool PetscStackActive(void) {return PETSC_FALSE;}
# 891 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscerror.h"
extern PetscErrorCode PetscStackCreate(void);
extern PetscErrorCode PetscStackView(FILE*);
extern PetscErrorCode PetscStackDestroy(void);
# 1560 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h" 2


extern PetscClassId PETSC_LARGEST_CLASSID;
extern PetscClassId PETSC_OBJECT_CLASSID;
extern PetscErrorCode PetscClassIdRegister(const char[],PetscClassId *);




extern PetscErrorCode PetscMemoryGetCurrentUsage(PetscLogDouble *);
extern PetscErrorCode PetscMemoryGetMaximumUsage(PetscLogDouble *);
extern PetscErrorCode PetscMemorySetGetMaximumUsage(void);
extern PetscErrorCode PetscMemoryTrace(const char[]);

extern PetscErrorCode PetscInfoAllow(PetscBool ,const char []);
extern PetscErrorCode PetscSleep(PetscReal);




extern PetscErrorCode PetscInitialize(int*,char***,const char[],const char[]);
extern PetscErrorCode PetscInitializeNoPointers(int,char**,const char[],const char[]);
extern PetscErrorCode PetscInitializeNoArguments(void);
extern PetscErrorCode PetscInitialized(PetscBool *);
extern PetscErrorCode PetscFinalized(PetscBool *);
extern PetscErrorCode PetscFinalize(void);
extern PetscErrorCode PetscInitializeFortran(void);
extern PetscErrorCode PetscGetArgs(int*,char ***);
extern PetscErrorCode PetscGetArguments(char ***);
extern PetscErrorCode PetscFreeArguments(char **);

extern PetscErrorCode PetscEnd(void);
extern PetscErrorCode PetscSysInitializePackage(void);

extern PetscErrorCode PetscPythonInitialize(const char[],const char[]);
extern PetscErrorCode PetscPythonFinalize(void);
extern PetscErrorCode PetscPythonPrintError(void);
extern PetscErrorCode PetscPythonMonitorSet(PetscObject,const char[]);





 typedef void (**PetscVoidStarFunction)(void);
 typedef void (*PetscVoidFunction)(void);
 typedef PetscErrorCode (*PetscErrorCodeFunction)(void);




extern PetscErrorCode PetscObjectDestroy(PetscObject*);
extern PetscErrorCode PetscObjectGetComm(PetscObject,MPI_Comm *);
extern PetscErrorCode PetscObjectGetClassId(PetscObject,PetscClassId *);
extern PetscErrorCode PetscObjectGetClassName(PetscObject,const char *[]);
extern PetscErrorCode PetscObjectSetType(PetscObject,const char []);
extern PetscErrorCode PetscObjectSetPrecision(PetscObject,PetscPrecision);
extern PetscErrorCode PetscObjectGetType(PetscObject,const char *[]);
extern PetscErrorCode PetscObjectSetName(PetscObject,const char[]);
extern PetscErrorCode PetscObjectGetName(PetscObject,const char*[]);
extern PetscErrorCode PetscObjectSetTabLevel(PetscObject,PetscInt);
extern PetscErrorCode PetscObjectGetTabLevel(PetscObject,PetscInt*);
extern PetscErrorCode PetscObjectIncrementTabLevel(PetscObject,PetscObject,PetscInt);
extern PetscErrorCode PetscObjectReference(PetscObject);
extern PetscErrorCode PetscObjectGetReference(PetscObject,PetscInt*);
extern PetscErrorCode PetscObjectDereference(PetscObject);
extern PetscErrorCode PetscObjectGetNewTag(PetscObject,PetscMPIInt *);
extern PetscErrorCode PetscObjectCompose(PetscObject,const char[],PetscObject);
extern PetscErrorCode PetscObjectRemoveReference(PetscObject,const char[]);
extern PetscErrorCode PetscObjectQuery(PetscObject,const char[],PetscObject *);
extern PetscErrorCode PetscObjectComposeFunction_Private(PetscObject,const char[],void (*)(void));

extern PetscErrorCode PetscObjectSetFromOptions(PetscObject);
extern PetscErrorCode PetscObjectSetUp(PetscObject);
extern PetscErrorCode PetscCommGetNewTag(MPI_Comm,PetscMPIInt *);
extern PetscErrorCode PetscObjectAddOptionsHandler(PetscObject,PetscErrorCode (*)(PetscObject,void*),PetscErrorCode (*)(PetscObject,void*),void*);
extern PetscErrorCode PetscObjectProcessOptionsHandlers(PetscObject);
extern PetscErrorCode PetscObjectDestroyOptionsHandlers(PetscObject);
extern PetscErrorCode PetscObjectsListGetGlobalNumbering(MPI_Comm,PetscInt,PetscObject*,PetscInt*,PetscInt*);

# 1 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscviewertypes.h" 1
# 18 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscviewertypes.h"
typedef struct _p_PetscViewer* PetscViewer;
# 1640 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h" 2
# 1 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscoptions.h" 1





# 1 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h" 1
# 7 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscoptions.h" 2


extern PetscErrorCode PetscOptionsHasName(const char[],const char[],PetscBool *);
extern PetscErrorCode PetscOptionsGetInt(const char[],const char [],PetscInt *,PetscBool *);
extern PetscErrorCode PetscOptionsGetBool(const char[],const char [],PetscBool *,PetscBool *);
extern PetscErrorCode PetscOptionsGetReal(const char[],const char[],PetscReal *,PetscBool *);
extern PetscErrorCode PetscOptionsGetScalar(const char[],const char[],PetscScalar *,PetscBool *);
extern PetscErrorCode PetscOptionsGetIntArray(const char[],const char[],PetscInt[],PetscInt *,PetscBool *);
extern PetscErrorCode PetscOptionsGetRealArray(const char[],const char[],PetscReal[],PetscInt *,PetscBool *);
extern PetscErrorCode PetscOptionsGetScalarArray(const char[],const char[],PetscScalar[],PetscInt *,PetscBool *);
extern PetscErrorCode PetscOptionsGetBoolArray(const char[],const char[],PetscBool [],PetscInt *,PetscBool *);
extern PetscErrorCode PetscOptionsGetString(const char[],const char[],char[],size_t,PetscBool *);
extern PetscErrorCode PetscOptionsGetStringArray(const char[],const char[],char*[],PetscInt*,PetscBool *);
extern PetscErrorCode PetscOptionsGetEList(const char[],const char[],const char*const*,PetscInt,PetscInt*,PetscBool *);
extern PetscErrorCode PetscOptionsGetEnum(const char[],const char[],const char*const*,PetscEnum*,PetscBool *);
extern PetscErrorCode PetscOptionsGetEnumArray(const char[],const char[],const char*const*,PetscEnum*,PetscInt *,PetscBool *);
extern PetscErrorCode PetscOptionsValidKey(const char[],PetscBool *);

extern PetscErrorCode PetscOptionsSetAlias(const char[],const char[]);
extern PetscErrorCode PetscOptionsSetValue(const char[],const char[]);
extern PetscErrorCode PetscOptionsClearValue(const char[]);

extern PetscErrorCode PetscOptionsAllUsed(PetscInt*);
extern PetscErrorCode PetscOptionsUsed(const char *,PetscBool*);
extern PetscErrorCode PetscOptionsLeft(void);
extern PetscErrorCode PetscOptionsView(PetscViewer);

extern PetscErrorCode PetscOptionsCreate(void);
extern PetscErrorCode PetscOptionsInsert(int*,char ***,const char[]);
extern PetscErrorCode PetscOptionsInsertFile(MPI_Comm,const char[],PetscBool );



extern PetscErrorCode PetscOptionsInsertString(const char[]);
extern PetscErrorCode PetscOptionsDestroy(void);
extern PetscErrorCode PetscOptionsClear(void);
extern PetscErrorCode PetscOptionsPrefixPush(const char[]);
extern PetscErrorCode PetscOptionsPrefixPop(void);

extern PetscErrorCode PetscOptionsReject(const char[],const char[]);
extern PetscErrorCode PetscOptionsGetAll(char*[]);

extern PetscErrorCode PetscOptionsGetenv(MPI_Comm,const char[],char[],size_t,PetscBool *);
extern PetscErrorCode PetscOptionsStringToInt(const char[],PetscInt*);
extern PetscErrorCode PetscOptionsStringToReal(const char[],PetscReal*);
extern PetscErrorCode PetscOptionsStringToBool(const char[],PetscBool*);

extern PetscErrorCode PetscOptionsMonitorSet(PetscErrorCode (*)(const char[], const char[], void*), void *, PetscErrorCode (*)(void**));
extern PetscErrorCode PetscOptionsMonitorCancel(void);
extern PetscErrorCode PetscOptionsMonitorDefault(const char[], const char[], void *);

extern PetscBool PetscOptionsPublish;





typedef enum {OPTION_INT,OPTION_BOOL,OPTION_REAL,OPTION_FLIST,OPTION_STRING,OPTION_REAL_ARRAY,OPTION_SCALAR_ARRAY,OPTION_HEAD,OPTION_INT_ARRAY,OPTION_ELIST,OPTION_BOOL_ARRAY,OPTION_STRING_ARRAY} PetscOptionType;
typedef struct _n_PetscOption* PetscOption;
struct _n_PetscOption{
  char *option;
  char *text;
  void *data;
  PetscFunctionList flist;
  const char *const *list;
  char nlist;
  char *man;
  size_t arraylength;
  PetscBool set;
  PetscOptionType type;
  PetscOption next;
  char *pman;
  void *edata;
};

typedef struct _p_PetscOptions {
  PetscInt count;
  PetscOption next;
  char *prefix,*pprefix;
  char *title;
  MPI_Comm comm;
  PetscBool printhelp,changedmethod,alreadyprinted;
  PetscObject object;
} PetscOptions;
# 202 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscoptions.h"
extern PetscErrorCode PetscOptionsBegin_Private(PetscOptions *,MPI_Comm,const char[],const char[],const char[]);
extern PetscErrorCode PetscObjectOptionsBegin_Private(PetscOptions *,PetscObject);
extern PetscErrorCode PetscOptionsEnd_Private(PetscOptions *);
extern PetscErrorCode PetscOptionsHead(PetscOptions *,const char[]);
# 260 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscoptions.h"
extern PetscErrorCode PetscOptionsEnum_Private(PetscOptions *,const char[],const char[],const char[],const char *const*,PetscEnum,PetscEnum*,PetscBool *);
extern PetscErrorCode PetscOptionsInt_Private(PetscOptions *,const char[],const char[],const char[],PetscInt,PetscInt*,PetscBool *);
extern PetscErrorCode PetscOptionsReal_Private(PetscOptions *,const char[],const char[],const char[],PetscReal,PetscReal*,PetscBool *);
extern PetscErrorCode PetscOptionsScalar_Private(PetscOptions *,const char[],const char[],const char[],PetscScalar,PetscScalar*,PetscBool *);
extern PetscErrorCode PetscOptionsName_Private(PetscOptions *,const char[],const char[],const char[],PetscBool *);
extern PetscErrorCode PetscOptionsString_Private(PetscOptions *,const char[],const char[],const char[],const char[],char*,size_t,PetscBool *);
extern PetscErrorCode PetscOptionsBool_Private(PetscOptions *,const char[],const char[],const char[],PetscBool ,PetscBool *,PetscBool *);
extern PetscErrorCode PetscOptionsBoolGroupBegin_Private(PetscOptions *,const char[],const char[],const char[],PetscBool *);
extern PetscErrorCode PetscOptionsBoolGroup_Private(PetscOptions *,const char[],const char[],const char[],PetscBool *);
extern PetscErrorCode PetscOptionsBoolGroupEnd_Private(PetscOptions *,const char[],const char[],const char[],PetscBool *);
extern PetscErrorCode PetscOptionsFList_Private(PetscOptions *,const char[],const char[],const char[],PetscFunctionList,const char[],char[],size_t,PetscBool *);
extern PetscErrorCode PetscOptionsEList_Private(PetscOptions *,const char[],const char[],const char[],const char*const*,PetscInt,const char[],PetscInt*,PetscBool *);
extern PetscErrorCode PetscOptionsRealArray_Private(PetscOptions *,const char[],const char[],const char[],PetscReal[],PetscInt*,PetscBool *);
extern PetscErrorCode PetscOptionsScalarArray_Private(PetscOptions *,const char[],const char[],const char[],PetscScalar[],PetscInt*,PetscBool *);
extern PetscErrorCode PetscOptionsIntArray_Private(PetscOptions *,const char[],const char[],const char[],PetscInt[],PetscInt*,PetscBool *);
extern PetscErrorCode PetscOptionsStringArray_Private(PetscOptions *,const char[],const char[],const char[],char*[],PetscInt*,PetscBool *);
extern PetscErrorCode PetscOptionsBoolArray_Private(PetscOptions *,const char[],const char[],const char[],PetscBool [],PetscInt*,PetscBool *);
extern PetscErrorCode PetscOptionsEnumArray_Private(PetscOptions *,const char[],const char[],const char[],const char *const*,PetscEnum[],PetscInt*,PetscBool *);


extern PetscErrorCode PetscOptionsSetFromOptions(void);
extern PetscErrorCode PetscOptionsSAWsDestroy(void);
# 1641 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h" 2

extern PetscErrorCode PetscMemoryShowUsage(PetscViewer,const char[]);
extern PetscErrorCode PetscObjectPrintClassNamePrefixType(PetscObject,PetscViewer);
extern PetscErrorCode PetscObjectView(PetscObject,PetscViewer);

extern PetscErrorCode PetscObjectQueryFunction_Private(PetscObject,const char[],void (**)(void));
extern PetscErrorCode PetscObjectSetOptionsPrefix(PetscObject,const char[]);
extern PetscErrorCode PetscObjectAppendOptionsPrefix(PetscObject,const char[]);
extern PetscErrorCode PetscObjectPrependOptionsPrefix(PetscObject,const char[]);
extern PetscErrorCode PetscObjectGetOptionsPrefix(PetscObject,const char*[]);
extern PetscErrorCode PetscObjectChangeTypeName(PetscObject,const char[]);
extern PetscErrorCode PetscObjectRegisterDestroy(PetscObject);
extern PetscErrorCode PetscObjectRegisterDestroyAll(void);
extern PetscErrorCode PetscObjectViewFromOptions(PetscObject,PetscObject,const char[]);
extern PetscErrorCode PetscObjectName(PetscObject);
extern PetscErrorCode PetscObjectTypeCompare(PetscObject,const char[],PetscBool *);
extern PetscErrorCode PetscObjectTypeCompareAny(PetscObject,PetscBool*,const char[],...);
extern PetscErrorCode PetscRegisterFinalize(PetscErrorCode (*)(void));
extern PetscErrorCode PetscRegisterFinalizeAll(void);
# 1687 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
typedef void* PetscDLHandle;
typedef enum {PETSC_DL_DECIDE=0,PETSC_DL_NOW=1,PETSC_DL_LOCAL=2} PetscDLMode;
extern PetscErrorCode PetscDLOpen(const char[],PetscDLMode,PetscDLHandle *);
extern PetscErrorCode PetscDLClose(PetscDLHandle *);
extern PetscErrorCode PetscDLSym(PetscDLHandle,const char[],void **);





extern PetscErrorCode PetscObjectsDump(FILE*,PetscBool);
# 1708 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
typedef struct _n_PetscObjectList *PetscObjectList;

extern PetscErrorCode PetscObjectListDestroy(PetscObjectList*);
extern PetscErrorCode PetscObjectListFind(PetscObjectList,const char[],PetscObject*);
extern PetscErrorCode PetscObjectListReverseFind(PetscObjectList,PetscObject,char**,PetscBool*);
extern PetscErrorCode PetscObjectListAdd(PetscObjectList *,const char[],PetscObject);
extern PetscErrorCode PetscObjectListRemoveReference(PetscObjectList *,const char[]);
extern PetscErrorCode PetscObjectListDuplicate(PetscObjectList,PetscObjectList *);







extern PetscErrorCode PetscFunctionListAdd_Private(PetscFunctionList*,const char[],void (*)(void));
extern PetscErrorCode PetscFunctionListDestroy(PetscFunctionList*);

extern PetscErrorCode PetscFunctionListFind_Private(PetscFunctionList,const char[],void (**)(void));
extern PetscErrorCode PetscFunctionListPrintTypes(MPI_Comm,FILE*,const char[],const char[],const char[],const char[],PetscFunctionList,const char[]);
extern PetscErrorCode PetscFunctionListDuplicate(PetscFunctionList,PetscFunctionList *);
extern PetscErrorCode PetscFunctionListView(PetscFunctionList,PetscViewer);
extern PetscErrorCode PetscFunctionListGet(PetscFunctionList,const char ***,int*);
# 1739 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
typedef struct _n_PetscDLLibrary *PetscDLLibrary;
extern PetscDLLibrary PetscDLLibrariesLoaded;
extern PetscErrorCode PetscDLLibraryAppend(MPI_Comm,PetscDLLibrary *,const char[]);
extern PetscErrorCode PetscDLLibraryPrepend(MPI_Comm,PetscDLLibrary *,const char[]);
extern PetscErrorCode PetscDLLibrarySym(MPI_Comm,PetscDLLibrary *,const char[],const char[],void **);
extern PetscErrorCode PetscDLLibraryPrintPath(PetscDLLibrary);
extern PetscErrorCode PetscDLLibraryRetrieve(MPI_Comm,const char[],char *,size_t,PetscBool *);
extern PetscErrorCode PetscDLLibraryOpen(MPI_Comm,const char[],PetscDLLibrary *);
extern PetscErrorCode PetscDLLibraryClose(PetscDLLibrary);




extern PetscErrorCode PetscSplitOwnership(MPI_Comm,PetscInt*,PetscInt*);
extern PetscErrorCode PetscSplitOwnershipBlock(MPI_Comm,PetscInt,PetscInt*,PetscInt*);
extern PetscErrorCode PetscSequentialPhaseBegin(MPI_Comm,PetscMPIInt);
extern PetscErrorCode PetscSequentialPhaseEnd(MPI_Comm,PetscMPIInt);
extern PetscErrorCode PetscBarrier(PetscObject);
extern PetscErrorCode PetscMPIDump(FILE*);
# 1791 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
extern PetscErrorCode (*PetscHelpPrintf)(MPI_Comm,const char[],...);




# 1 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petsclog.h" 1
# 17 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petsclog.h"
typedef int PetscLogEvent;
# 26 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petsclog.h"
typedef int PetscLogStage;


extern PetscLogEvent PETSC_LARGEST_EVENT;


extern PetscLogDouble petsc_TotalFlops;
extern PetscLogDouble petsc_tmp_flops;


extern PetscErrorCode PetscInfo_Private(const char[],void*,const char[],...);
# 56 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petsclog.h"
extern PetscErrorCode PetscInfoDeactivateClass(PetscClassId);
extern PetscErrorCode PetscInfoActivateClass(PetscClassId);
extern PetscBool PetscLogPrintInfo;







typedef struct _n_PetscIntStack *PetscIntStack;
# 75 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petsclog.h"
typedef struct {
  char *name;
  PetscClassId classid;
} PetscClassRegInfo;

typedef struct {
  PetscClassId id;
  int creations;
  int destructions;
  PetscLogDouble mem;
  PetscLogDouble descMem;
} PetscClassPerfInfo;

typedef struct _n_PetscClassRegLog *PetscClassRegLog;
struct _n_PetscClassRegLog {
  int numClasses;
  int maxClasses;
  PetscClassRegInfo *classInfo;
};

typedef struct _n_PetscClassPerfLog *PetscClassPerfLog;
struct _n_PetscClassPerfLog {
  int numClasses;
  int maxClasses;
  PetscClassPerfInfo *classInfo;
};
# 111 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petsclog.h"
typedef struct {
  char *name;
  PetscClassId classid;




} PetscEventRegInfo;

typedef struct {
  int id;
  PetscBool active;
  PetscBool visible;
  int depth;
  int count;
  PetscLogDouble flops, flops2,flopsTmp;
  PetscLogDouble time, time2, timeTmp;
  PetscLogDouble numMessages;
  PetscLogDouble messageLength;
  PetscLogDouble numReductions;
} PetscEventPerfInfo;

typedef struct _n_PetscEventRegLog *PetscEventRegLog;
struct _n_PetscEventRegLog {
  int numEvents;
  int maxEvents;
  PetscEventRegInfo *eventInfo;
};

typedef struct _n_PetscEventPerfLog *PetscEventPerfLog;
struct _n_PetscEventPerfLog {
  int numEvents;
  int maxEvents;
  PetscEventPerfInfo *eventInfo;
};






typedef struct _PetscStageInfo {
  char *name;
  PetscBool used;
  PetscEventPerfInfo perfInfo;
  PetscEventPerfLog eventLog;
  PetscClassPerfLog classLog;
} PetscStageInfo;

typedef struct _n_PetscStageLog *PetscStageLog;
struct _n_PetscStageLog {
  int numStages;
  int maxStages;
  PetscIntStack stack;
  int curStage;
  PetscStageInfo *stageInfo;
  PetscEventRegLog eventLog;
  PetscClassRegLog classLog;
};

extern PetscErrorCode PetscLogGetStageLog(PetscStageLog*);
extern PetscErrorCode PetscStageLogGetCurrent(PetscStageLog,int*);
extern PetscErrorCode PetscStageLogGetEventPerfLog(PetscStageLog,int,PetscEventPerfLog*);

extern PetscErrorCode PetscLogObjectParent(PetscObject,PetscObject);
extern PetscErrorCode PetscLogObjectMemory(PetscObject,PetscLogDouble);



extern PetscStageLog petsc_stageLog;
# 203 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petsclog.h"
static inline PetscErrorCode PetscLogFlops(PetscLogDouble n)
{
  ;



  petsc_TotalFlops += 1.0*n;
  return(0);
}






extern PetscErrorCode (*PetscLogPLB)(PetscLogEvent,int,PetscObject,PetscObject,PetscObject,PetscObject);
extern PetscErrorCode (*PetscLogPLE)(PetscLogEvent,int,PetscObject,PetscObject,PetscObject,PetscObject);
extern PetscErrorCode (*PetscLogPHC)(PetscObject);
extern PetscErrorCode (*PetscLogPHD)(PetscObject);





extern PetscErrorCode PetscLogBegin(void);
extern PetscErrorCode PetscLogAllBegin(void);
extern PetscErrorCode PetscLogTraceBegin(FILE *);
extern PetscErrorCode PetscLogActions(PetscBool);
extern PetscErrorCode PetscLogObjects(PetscBool);

extern PetscErrorCode PetscLogDestroy(void);
extern PetscErrorCode PetscLogSet(PetscErrorCode (*)(int, int, PetscObject, PetscObject, PetscObject, PetscObject),
                                   PetscErrorCode (*)(int, int, PetscObject, PetscObject, PetscObject, PetscObject));
extern PetscErrorCode PetscLogObjectState(PetscObject, const char[], ...);

extern PetscErrorCode PetscLogView(PetscViewer);
extern PetscErrorCode PetscLogViewFromOptions(void);
extern PetscErrorCode PetscLogDump(const char[]);

extern PetscErrorCode PetscGetFlops(PetscLogDouble *);

extern PetscErrorCode PetscLogStageRegister(const char[],PetscLogStage*);
extern PetscErrorCode PetscLogStagePush(PetscLogStage);
extern PetscErrorCode PetscLogStagePop(void);
extern PetscErrorCode PetscLogStageSetActive(PetscLogStage, PetscBool );
extern PetscErrorCode PetscLogStageGetActive(PetscLogStage, PetscBool *);
extern PetscErrorCode PetscLogStageSetVisible(PetscLogStage, PetscBool );
extern PetscErrorCode PetscLogStageGetVisible(PetscLogStage, PetscBool *);
extern PetscErrorCode PetscLogStageGetId(const char [], PetscLogStage *);

extern PetscErrorCode PetscLogEventRegister(const char[], PetscClassId,PetscLogEvent*);
extern PetscErrorCode PetscLogEventActivate(PetscLogEvent);
extern PetscErrorCode PetscLogEventDeactivate(PetscLogEvent);
extern PetscErrorCode PetscLogEventSetActiveAll(PetscLogEvent, PetscBool );
extern PetscErrorCode PetscLogEventActivateClass(PetscClassId);
extern PetscErrorCode PetscLogEventDeactivateClass(PetscClassId);
extern PetscErrorCode PetscLogEventGetId(const char[],PetscLogEvent*);
extern PetscErrorCode PetscLogEventGetPerfInfo(int, PetscLogEvent, PetscEventPerfInfo *);


extern PetscLogDouble petsc_irecv_ct;
extern PetscLogDouble petsc_isend_ct;
extern PetscLogDouble petsc_recv_ct;
extern PetscLogDouble petsc_send_ct;
extern PetscLogDouble petsc_irecv_len;
extern PetscLogDouble petsc_isend_len;
extern PetscLogDouble petsc_recv_len;
extern PetscLogDouble petsc_send_len;
extern PetscLogDouble petsc_allreduce_ct;
extern PetscLogDouble petsc_gather_ct;
extern PetscLogDouble petsc_scatter_ct;
extern PetscLogDouble petsc_wait_ct;
extern PetscLogDouble petsc_wait_any_ct;
extern PetscLogDouble petsc_wait_all_ct;
extern PetscLogDouble petsc_sum_of_waits_ct;
# 294 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petsclog.h"
extern PetscErrorCode PetscLogEventGetFlops(PetscLogEvent, PetscLogDouble*);
extern PetscErrorCode PetscLogEventZeroFlops(PetscLogEvent);
# 314 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petsclog.h"
static inline PetscErrorCode PetscMPITypeSize(PetscLogDouble *buff,PetscMPIInt count,MPI_Datatype type)
{
  PetscMPIInt mysize;
  if (type == ((MPI_Datatype)0x0c000000)) return 0;
  else return (MPI_Type_size(type,&mysize) || ((*buff += (PetscLogDouble) (count*mysize)),0));
}

static inline PetscErrorCode PetscMPITypeSizeComm(MPI_Comm comm, PetscLogDouble *buff,PetscMPIInt *counts,MPI_Datatype type)
{
  PetscMPIInt mysize, commsize, p;
  PetscErrorCode _myierr;

  if (type == ((MPI_Datatype)0x0c000000)) return 0;
  _myierr = MPI_Comm_size(comm,&commsize);do {if (__builtin_expect(!!(_myierr),0)) return PetscError(((MPI_Comm)0x44000001),327,__func__,"/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petsclog.h",_myierr,PETSC_ERROR_REPEAT," ");} while (0);
  _myierr = MPI_Type_size(type,&mysize);do {if (__builtin_expect(!!(_myierr),0)) return PetscError(((MPI_Comm)0x44000001),328,__func__,"/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petsclog.h",_myierr,PETSC_ERROR_REPEAT," ");} while (0);
  for (p = 0; p < commsize; ++p) {
    *buff += (PetscLogDouble) (counts[p]*mysize);
  }
  return 0;
}




static inline int PetscMPIParallelComm(MPI_Comm comm)
{
  PetscMPIInt size; MPI_Comm_size(comm,&size); return size > 1;
}
# 465 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petsclog.h"
extern PetscErrorCode PetscIntStackCreate(PetscIntStack *);
extern PetscErrorCode PetscIntStackDestroy(PetscIntStack);
extern PetscErrorCode PetscIntStackPush(PetscIntStack, int);
extern PetscErrorCode PetscIntStackPop(PetscIntStack, int *);
extern PetscErrorCode PetscIntStackTop(PetscIntStack, int *);
extern PetscErrorCode PetscIntStackEmpty(PetscIntStack, PetscBool *);
# 510 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petsclog.h"
extern PetscBool PetscPreLoadingUsed;
extern PetscBool PetscPreLoadingOn;
# 1797 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h" 2






extern PetscErrorCode PetscFixFilename(const char[],char[]);
extern PetscErrorCode PetscFOpen(MPI_Comm,const char[],const char[],FILE**);
extern PetscErrorCode PetscFClose(MPI_Comm,FILE*);
extern PetscErrorCode PetscFPrintf(MPI_Comm,FILE*,const char[],...);
extern PetscErrorCode PetscPrintf(MPI_Comm,const char[],...);
extern PetscErrorCode PetscSNPrintf(char*,size_t,const char [],...);
extern PetscErrorCode PetscSNPrintfCount(char*,size_t,const char [],size_t*,...);


# 1 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/stdarg.h" 1 3 4
# 1813 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h" 2
extern PetscErrorCode PetscVSNPrintf(char*,size_t,const char[],size_t*,va_list);
extern PetscErrorCode (*PetscVFPrintf)(FILE*,const char[],va_list);
extern PetscErrorCode PetscVFPrintfDefault(FILE*,const char[],va_list);
# 1825 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
extern PetscErrorCode PetscErrorPrintfDefault(const char [],...);
extern PetscErrorCode PetscErrorPrintfNone(const char [],...);
extern PetscErrorCode PetscHelpPrintfDefault(MPI_Comm,const char [],...);


extern PetscErrorCode PetscPOpen(MPI_Comm,const char[],const char[],const char[],FILE **);
extern PetscErrorCode PetscPClose(MPI_Comm,FILE*,int*);
extern PetscErrorCode PetscPOpenSetMachine(const char[]);


extern PetscErrorCode PetscSynchronizedPrintf(MPI_Comm,const char[],...);
extern PetscErrorCode PetscSynchronizedFPrintf(MPI_Comm,FILE*,const char[],...);
extern PetscErrorCode PetscSynchronizedFlush(MPI_Comm,FILE*);
extern PetscErrorCode PetscSynchronizedFGets(MPI_Comm,FILE*,size_t,char[]);
extern PetscErrorCode PetscStartMatlab(MPI_Comm,const char[],const char[],FILE**);
extern PetscErrorCode PetscStartJava(MPI_Comm,const char[],const char[],FILE**);
extern PetscErrorCode PetscGetPetscDir(const char*[]);

extern PetscErrorCode PetscPopUpSelect(MPI_Comm,const char*,const char*,int,const char**,int*);
# 1852 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
extern PetscClassId PETSC_CONTAINER_CLASSID;
typedef struct _p_PetscContainer* PetscContainer;
extern PetscErrorCode PetscContainerGetPointer(PetscContainer,void **);
extern PetscErrorCode PetscContainerSetPointer(PetscContainer,void *);
extern PetscErrorCode PetscContainerDestroy(PetscContainer*);
extern PetscErrorCode PetscContainerCreate(MPI_Comm,PetscContainer *);
extern PetscErrorCode PetscContainerSetUserDestroy(PetscContainer, PetscErrorCode (*)(void*));




extern PetscMPIInt PetscGlobalRank;
extern PetscMPIInt PetscGlobalSize;
extern PetscErrorCode PetscIntView(PetscInt,const PetscInt[],PetscViewer);
extern PetscErrorCode PetscRealView(PetscInt,const PetscReal[],PetscViewer);
extern PetscErrorCode PetscScalarView(PetscInt,const PetscScalar[],PetscViewer);

# 1 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/stddef.h" 1 3 4
# 149 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/stddef.h" 3 4

# 149 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/stddef.h" 3 4
typedef long int ptrdiff_t;
# 328 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/stddef.h" 3 4
typedef int wchar_t;
# 426 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/stddef.h" 3 4
typedef struct {
  long long __max_align_ll __attribute__((__aligned__(__alignof__(long long))));
  long double __max_align_ld __attribute__((__aligned__(__alignof__(long double))));
# 437 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/stddef.h" 3 4
} max_align_t;
# 1870 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h" 2
# 1 "/usr/include/string.h" 1 3 4
# 27 "/usr/include/string.h" 3 4





# 1 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/stddef.h" 1 3 4
# 33 "/usr/include/string.h" 2 3 4









extern void *memcpy (void *__restrict __dest, const void *__restrict __src,
       size_t __n) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));


extern void *memmove (void *__dest, const void *__src, size_t __n)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));






extern void *memccpy (void *__restrict __dest, const void *__restrict __src,
        int __c, size_t __n)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));





extern void *memset (void *__s, int __c, size_t __n) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));


extern int memcmp (const void *__s1, const void *__s2, size_t __n)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1, 2)));
# 92 "/usr/include/string.h" 3 4
extern void *memchr (const void *__s, int __c, size_t __n)
      __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1)));


# 123 "/usr/include/string.h" 3 4


extern char *strcpy (char *__restrict __dest, const char *__restrict __src)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));

extern char *strncpy (char *__restrict __dest,
        const char *__restrict __src, size_t __n)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));


extern char *strcat (char *__restrict __dest, const char *__restrict __src)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));

extern char *strncat (char *__restrict __dest, const char *__restrict __src,
        size_t __n) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));


extern int strcmp (const char *__s1, const char *__s2)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1, 2)));

extern int strncmp (const char *__s1, const char *__s2, size_t __n)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1, 2)));


extern int strcoll (const char *__s1, const char *__s2)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1, 2)));

extern size_t strxfrm (char *__restrict __dest,
         const char *__restrict __src, size_t __n)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (2)));






# 1 "/usr/include/xlocale.h" 1 3 4
# 27 "/usr/include/xlocale.h" 3 4
typedef struct __locale_struct
{

  struct __locale_data *__locales[13];


  const unsigned short int *__ctype_b;
  const int *__ctype_tolower;
  const int *__ctype_toupper;


  const char *__names[13];
} *__locale_t;


typedef __locale_t locale_t;
# 160 "/usr/include/string.h" 2 3 4


extern int strcoll_l (const char *__s1, const char *__s2, __locale_t __l)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1, 2, 3)));

extern size_t strxfrm_l (char *__dest, const char *__src, size_t __n,
    __locale_t __l) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (2, 4)));





extern char *strdup (const char *__s)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__malloc__)) __attribute__ ((__nonnull__ (1)));






extern char *strndup (const char *__string, size_t __n)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__malloc__)) __attribute__ ((__nonnull__ (1)));
# 207 "/usr/include/string.h" 3 4

# 232 "/usr/include/string.h" 3 4
extern char *strchr (const char *__s, int __c)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1)));
# 259 "/usr/include/string.h" 3 4
extern char *strrchr (const char *__s, int __c)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1)));


# 278 "/usr/include/string.h" 3 4



extern size_t strcspn (const char *__s, const char *__reject)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1, 2)));


extern size_t strspn (const char *__s, const char *__accept)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1, 2)));
# 311 "/usr/include/string.h" 3 4
extern char *strpbrk (const char *__s, const char *__accept)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1, 2)));
# 338 "/usr/include/string.h" 3 4
extern char *strstr (const char *__haystack, const char *__needle)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1, 2)));




extern char *strtok (char *__restrict __s, const char *__restrict __delim)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (2)));




extern char *__strtok_r (char *__restrict __s,
    const char *__restrict __delim,
    char **__restrict __save_ptr)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (2, 3)));

extern char *strtok_r (char *__restrict __s, const char *__restrict __delim,
         char **__restrict __save_ptr)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (2, 3)));
# 393 "/usr/include/string.h" 3 4


extern size_t strlen (const char *__s)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1)));





extern size_t strnlen (const char *__string, size_t __maxlen)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1)));





extern char *strerror (int __errnum) __attribute__ ((__nothrow__ , __leaf__));

# 423 "/usr/include/string.h" 3 4
extern int strerror_r (int __errnum, char *__buf, size_t __buflen) __asm__ ("" "__xpg_strerror_r") __attribute__ ((__nothrow__ , __leaf__))

                        __attribute__ ((__nonnull__ (2)));
# 441 "/usr/include/string.h" 3 4
extern char *strerror_l (int __errnum, __locale_t __l) __attribute__ ((__nothrow__ , __leaf__));





extern void __bzero (void *__s, size_t __n) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));



extern void bcopy (const void *__src, void *__dest, size_t __n)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));


extern void bzero (void *__s, size_t __n) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));


extern int bcmp (const void *__s1, const void *__s2, size_t __n)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1, 2)));
# 485 "/usr/include/string.h" 3 4
extern char *index (const char *__s, int __c)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1)));
# 513 "/usr/include/string.h" 3 4
extern char *rindex (const char *__s, int __c)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1)));




extern int ffs (int __i) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
# 532 "/usr/include/string.h" 3 4
extern int strcasecmp (const char *__s1, const char *__s2)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1, 2)));


extern int strncasecmp (const char *__s1, const char *__s2, size_t __n)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1, 2)));
# 555 "/usr/include/string.h" 3 4
extern char *strsep (char **__restrict __stringp,
       const char *__restrict __delim)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));




extern char *strsignal (int __sig) __attribute__ ((__nothrow__ , __leaf__));


extern char *__stpcpy (char *__restrict __dest, const char *__restrict __src)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));
extern char *stpcpy (char *__restrict __dest, const char *__restrict __src)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));



extern char *__stpncpy (char *__restrict __dest,
   const char *__restrict __src, size_t __n)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));
extern char *stpncpy (char *__restrict __dest,
        const char *__restrict __src, size_t __n)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));
# 642 "/usr/include/string.h" 3 4

# 1871 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h" 2

# 1 "/usr/include/stdlib.h" 1 3 4
# 32 "/usr/include/stdlib.h" 3 4
# 1 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/stddef.h" 1 3 4
# 33 "/usr/include/stdlib.h" 2 3 4








# 1 "/usr/include/bits/waitflags.h" 1 3 4
# 42 "/usr/include/stdlib.h" 2 3 4
# 1 "/usr/include/bits/waitstatus.h" 1 3 4
# 64 "/usr/include/bits/waitstatus.h" 3 4
# 1 "/usr/include/endian.h" 1 3 4
# 36 "/usr/include/endian.h" 3 4
# 1 "/usr/include/bits/endian.h" 1 3 4
# 37 "/usr/include/endian.h" 2 3 4
# 60 "/usr/include/endian.h" 3 4
# 1 "/usr/include/bits/byteswap.h" 1 3 4
# 28 "/usr/include/bits/byteswap.h" 3 4
# 1 "/usr/include/bits/wordsize.h" 1 3 4
# 29 "/usr/include/bits/byteswap.h" 2 3 4






# 1 "/usr/include/bits/byteswap-16.h" 1 3 4
# 36 "/usr/include/bits/byteswap.h" 2 3 4
# 44 "/usr/include/bits/byteswap.h" 3 4
static __inline unsigned int
__bswap_32 (unsigned int __bsx)
{
  return __builtin_bswap32 (__bsx);
}
# 108 "/usr/include/bits/byteswap.h" 3 4
static __inline __uint64_t
__bswap_64 (__uint64_t __bsx)
{
  return __builtin_bswap64 (__bsx);
}
# 61 "/usr/include/endian.h" 2 3 4
# 65 "/usr/include/bits/waitstatus.h" 2 3 4

union wait
  {
    int w_status;
    struct
      {

 unsigned int __w_termsig:7;
 unsigned int __w_coredump:1;
 unsigned int __w_retcode:8;
 unsigned int:16;







      } __wait_terminated;
    struct
      {

 unsigned int __w_stopval:8;
 unsigned int __w_stopsig:8;
 unsigned int:16;






      } __wait_stopped;
  };
# 43 "/usr/include/stdlib.h" 2 3 4
# 67 "/usr/include/stdlib.h" 3 4
typedef union
  {
    union wait *__uptr;
    int *__iptr;
  } __WAIT_STATUS __attribute__ ((__transparent_union__));
# 95 "/usr/include/stdlib.h" 3 4


typedef struct
  {
    int quot;
    int rem;
  } div_t;



typedef struct
  {
    long int quot;
    long int rem;
  } ldiv_t;







__extension__ typedef struct
  {
    long long int quot;
    long long int rem;
  } lldiv_t;


# 139 "/usr/include/stdlib.h" 3 4
extern size_t __ctype_get_mb_cur_max (void) __attribute__ ((__nothrow__ , __leaf__)) ;




extern double atof (const char *__nptr)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1))) ;

extern int atoi (const char *__nptr)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1))) ;

extern long int atol (const char *__nptr)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1))) ;





__extension__ extern long long int atoll (const char *__nptr)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1))) ;





extern double strtod (const char *__restrict __nptr,
        char **__restrict __endptr)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));





extern float strtof (const char *__restrict __nptr,
       char **__restrict __endptr) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));

extern long double strtold (const char *__restrict __nptr,
       char **__restrict __endptr)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));





extern long int strtol (const char *__restrict __nptr,
   char **__restrict __endptr, int __base)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));

extern unsigned long int strtoul (const char *__restrict __nptr,
      char **__restrict __endptr, int __base)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));




__extension__
extern long long int strtoq (const char *__restrict __nptr,
        char **__restrict __endptr, int __base)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));

__extension__
extern unsigned long long int strtouq (const char *__restrict __nptr,
           char **__restrict __endptr, int __base)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));





__extension__
extern long long int strtoll (const char *__restrict __nptr,
         char **__restrict __endptr, int __base)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));

__extension__
extern unsigned long long int strtoull (const char *__restrict __nptr,
     char **__restrict __endptr, int __base)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));

# 305 "/usr/include/stdlib.h" 3 4
extern char *l64a (long int __n) __attribute__ ((__nothrow__ , __leaf__)) ;


extern long int a64l (const char *__s)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1))) ;




# 1 "/usr/include/sys/types.h" 1 3 4
# 27 "/usr/include/sys/types.h" 3 4






typedef __u_char u_char;
typedef __u_short u_short;
typedef __u_int u_int;
typedef __u_long u_long;
typedef __quad_t quad_t;
typedef __u_quad_t u_quad_t;
typedef __fsid_t fsid_t;




typedef __loff_t loff_t;



typedef __ino_t ino_t;
# 60 "/usr/include/sys/types.h" 3 4
typedef __dev_t dev_t;




typedef __gid_t gid_t;




typedef __mode_t mode_t;




typedef __nlink_t nlink_t;




typedef __uid_t uid_t;
# 98 "/usr/include/sys/types.h" 3 4
typedef __pid_t pid_t;





typedef __id_t id_t;
# 115 "/usr/include/sys/types.h" 3 4
typedef __daddr_t daddr_t;
typedef __caddr_t caddr_t;





typedef __key_t key_t;
# 132 "/usr/include/sys/types.h" 3 4
# 1 "/usr/include/time.h" 1 3 4
# 57 "/usr/include/time.h" 3 4


typedef __clock_t clock_t;



# 73 "/usr/include/time.h" 3 4


typedef __time_t time_t;



# 91 "/usr/include/time.h" 3 4
typedef __clockid_t clockid_t;
# 103 "/usr/include/time.h" 3 4
typedef __timer_t timer_t;
# 133 "/usr/include/sys/types.h" 2 3 4
# 146 "/usr/include/sys/types.h" 3 4
# 1 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/stddef.h" 1 3 4
# 147 "/usr/include/sys/types.h" 2 3 4



typedef unsigned long int ulong;
typedef unsigned short int ushort;
typedef unsigned int uint;
# 200 "/usr/include/sys/types.h" 3 4
typedef unsigned int u_int8_t __attribute__ ((__mode__ (__QI__)));
typedef unsigned int u_int16_t __attribute__ ((__mode__ (__HI__)));
typedef unsigned int u_int32_t __attribute__ ((__mode__ (__SI__)));
typedef unsigned int u_int64_t __attribute__ ((__mode__ (__DI__)));

typedef int register_t __attribute__ ((__mode__ (__word__)));
# 219 "/usr/include/sys/types.h" 3 4
# 1 "/usr/include/sys/select.h" 1 3 4
# 30 "/usr/include/sys/select.h" 3 4
# 1 "/usr/include/bits/select.h" 1 3 4
# 22 "/usr/include/bits/select.h" 3 4
# 1 "/usr/include/bits/wordsize.h" 1 3 4
# 23 "/usr/include/bits/select.h" 2 3 4
# 31 "/usr/include/sys/select.h" 2 3 4


# 1 "/usr/include/bits/sigset.h" 1 3 4
# 23 "/usr/include/bits/sigset.h" 3 4
typedef int __sig_atomic_t;




typedef struct
  {
    unsigned long int __val[(1024 / (8 * sizeof (unsigned long int)))];
  } __sigset_t;
# 34 "/usr/include/sys/select.h" 2 3 4



typedef __sigset_t sigset_t;





# 1 "/usr/include/time.h" 1 3 4
# 120 "/usr/include/time.h" 3 4
struct timespec
  {
    __time_t tv_sec;
    __syscall_slong_t tv_nsec;
  };
# 44 "/usr/include/sys/select.h" 2 3 4

# 1 "/usr/include/bits/time.h" 1 3 4
# 30 "/usr/include/bits/time.h" 3 4
struct timeval
  {
    __time_t tv_sec;
    __suseconds_t tv_usec;
  };
# 46 "/usr/include/sys/select.h" 2 3 4


typedef __suseconds_t suseconds_t;





typedef long int __fd_mask;
# 64 "/usr/include/sys/select.h" 3 4
typedef struct
  {






    __fd_mask __fds_bits[1024 / (8 * (int) sizeof (__fd_mask))];


  } fd_set;






typedef __fd_mask fd_mask;
# 96 "/usr/include/sys/select.h" 3 4

# 106 "/usr/include/sys/select.h" 3 4
extern int select (int __nfds, fd_set *__restrict __readfds,
     fd_set *__restrict __writefds,
     fd_set *__restrict __exceptfds,
     struct timeval *__restrict __timeout);
# 118 "/usr/include/sys/select.h" 3 4
extern int pselect (int __nfds, fd_set *__restrict __readfds,
      fd_set *__restrict __writefds,
      fd_set *__restrict __exceptfds,
      const struct timespec *__restrict __timeout,
      const __sigset_t *__restrict __sigmask);
# 131 "/usr/include/sys/select.h" 3 4

# 220 "/usr/include/sys/types.h" 2 3 4


# 1 "/usr/include/sys/sysmacros.h" 1 3 4
# 29 "/usr/include/sys/sysmacros.h" 3 4


__extension__
extern unsigned int gnu_dev_major (unsigned long long int __dev)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
__extension__
extern unsigned int gnu_dev_minor (unsigned long long int __dev)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
__extension__
extern unsigned long long int gnu_dev_makedev (unsigned int __major,
            unsigned int __minor)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
# 63 "/usr/include/sys/sysmacros.h" 3 4

# 223 "/usr/include/sys/types.h" 2 3 4





typedef __blksize_t blksize_t;






typedef __blkcnt_t blkcnt_t;



typedef __fsblkcnt_t fsblkcnt_t;



typedef __fsfilcnt_t fsfilcnt_t;
# 270 "/usr/include/sys/types.h" 3 4
# 1 "/usr/include/bits/pthreadtypes.h" 1 3 4
# 21 "/usr/include/bits/pthreadtypes.h" 3 4
# 1 "/usr/include/bits/wordsize.h" 1 3 4
# 22 "/usr/include/bits/pthreadtypes.h" 2 3 4
# 60 "/usr/include/bits/pthreadtypes.h" 3 4
typedef unsigned long int pthread_t;


union pthread_attr_t
{
  char __size[56];
  long int __align;
};

typedef union pthread_attr_t pthread_attr_t;





typedef struct __pthread_internal_list
{
  struct __pthread_internal_list *__prev;
  struct __pthread_internal_list *__next;
} __pthread_list_t;
# 90 "/usr/include/bits/pthreadtypes.h" 3 4
typedef union
{
  struct __pthread_mutex_s
  {
    int __lock;
    unsigned int __count;
    int __owner;

    unsigned int __nusers;



    int __kind;

    short __spins;
    short __elision;
    __pthread_list_t __list;
# 125 "/usr/include/bits/pthreadtypes.h" 3 4
  } __data;
  char __size[40];
  long int __align;
} pthread_mutex_t;

typedef union
{
  char __size[4];
  int __align;
} pthread_mutexattr_t;




typedef union
{
  struct
  {
    int __lock;
    unsigned int __futex;
    __extension__ unsigned long long int __total_seq;
    __extension__ unsigned long long int __wakeup_seq;
    __extension__ unsigned long long int __woken_seq;
    void *__mutex;
    unsigned int __nwaiters;
    unsigned int __broadcast_seq;
  } __data;
  char __size[48];
  __extension__ long long int __align;
} pthread_cond_t;

typedef union
{
  char __size[4];
  int __align;
} pthread_condattr_t;



typedef unsigned int pthread_key_t;



typedef int pthread_once_t;





typedef union
{

  struct
  {
    int __lock;
    unsigned int __nr_readers;
    unsigned int __readers_wakeup;
    unsigned int __writer_wakeup;
    unsigned int __nr_readers_queued;
    unsigned int __nr_writers_queued;
    int __writer;
    int __shared;
    unsigned long int __pad1;
    unsigned long int __pad2;


    unsigned int __flags;

  } __data;
# 212 "/usr/include/bits/pthreadtypes.h" 3 4
  char __size[56];
  long int __align;
} pthread_rwlock_t;

typedef union
{
  char __size[8];
  long int __align;
} pthread_rwlockattr_t;





typedef volatile int pthread_spinlock_t;




typedef union
{
  char __size[32];
  long int __align;
} pthread_barrier_t;

typedef union
{
  char __size[4];
  int __align;
} pthread_barrierattr_t;
# 271 "/usr/include/sys/types.h" 2 3 4



# 315 "/usr/include/stdlib.h" 2 3 4






extern long int random (void) __attribute__ ((__nothrow__ , __leaf__));


extern void srandom (unsigned int __seed) __attribute__ ((__nothrow__ , __leaf__));





extern char *initstate (unsigned int __seed, char *__statebuf,
   size_t __statelen) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (2)));



extern char *setstate (char *__statebuf) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));







struct random_data
  {
    int32_t *fptr;
    int32_t *rptr;
    int32_t *state;
    int rand_type;
    int rand_deg;
    int rand_sep;
    int32_t *end_ptr;
  };

extern int random_r (struct random_data *__restrict __buf,
       int32_t *__restrict __result) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));

extern int srandom_r (unsigned int __seed, struct random_data *__buf)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (2)));

extern int initstate_r (unsigned int __seed, char *__restrict __statebuf,
   size_t __statelen,
   struct random_data *__restrict __buf)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (2, 4)));

extern int setstate_r (char *__restrict __statebuf,
         struct random_data *__restrict __buf)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));






extern int rand (void) __attribute__ ((__nothrow__ , __leaf__));

extern void srand (unsigned int __seed) __attribute__ ((__nothrow__ , __leaf__));




extern int rand_r (unsigned int *__seed) __attribute__ ((__nothrow__ , __leaf__));







extern double drand48 (void) __attribute__ ((__nothrow__ , __leaf__));
extern double erand48 (unsigned short int __xsubi[3]) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));


extern long int lrand48 (void) __attribute__ ((__nothrow__ , __leaf__));
extern long int nrand48 (unsigned short int __xsubi[3])
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));


extern long int mrand48 (void) __attribute__ ((__nothrow__ , __leaf__));
extern long int jrand48 (unsigned short int __xsubi[3])
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));


extern void srand48 (long int __seedval) __attribute__ ((__nothrow__ , __leaf__));
extern unsigned short int *seed48 (unsigned short int __seed16v[3])
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));
extern void lcong48 (unsigned short int __param[7]) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));





struct drand48_data
  {
    unsigned short int __x[3];
    unsigned short int __old_x[3];
    unsigned short int __c;
    unsigned short int __init;
    unsigned long long int __a;
  };


extern int drand48_r (struct drand48_data *__restrict __buffer,
        double *__restrict __result) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));
extern int erand48_r (unsigned short int __xsubi[3],
        struct drand48_data *__restrict __buffer,
        double *__restrict __result) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));


extern int lrand48_r (struct drand48_data *__restrict __buffer,
        long int *__restrict __result)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));
extern int nrand48_r (unsigned short int __xsubi[3],
        struct drand48_data *__restrict __buffer,
        long int *__restrict __result)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));


extern int mrand48_r (struct drand48_data *__restrict __buffer,
        long int *__restrict __result)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));
extern int jrand48_r (unsigned short int __xsubi[3],
        struct drand48_data *__restrict __buffer,
        long int *__restrict __result)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));


extern int srand48_r (long int __seedval, struct drand48_data *__buffer)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (2)));

extern int seed48_r (unsigned short int __seed16v[3],
       struct drand48_data *__buffer) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));

extern int lcong48_r (unsigned short int __param[7],
        struct drand48_data *__buffer)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));









extern void *malloc (size_t __size) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__malloc__)) ;

extern void *calloc (size_t __nmemb, size_t __size)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__malloc__)) ;










extern void *realloc (void *__ptr, size_t __size)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__warn_unused_result__));

extern void free (void *__ptr) __attribute__ ((__nothrow__ , __leaf__));




extern void cfree (void *__ptr) __attribute__ ((__nothrow__ , __leaf__));



# 1 "/usr/include/alloca.h" 1 3 4
# 24 "/usr/include/alloca.h" 3 4
# 1 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/stddef.h" 1 3 4
# 25 "/usr/include/alloca.h" 2 3 4







extern void *alloca (size_t __size) __attribute__ ((__nothrow__ , __leaf__));






# 492 "/usr/include/stdlib.h" 2 3 4





extern void *valloc (size_t __size) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__malloc__)) ;




extern int posix_memalign (void **__memptr, size_t __alignment, size_t __size)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1))) ;




extern void *aligned_alloc (size_t __alignment, size_t __size)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__malloc__, __alloc_size__ (2)));




extern void abort (void) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__noreturn__));



extern int atexit (void (*__func) (void)) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));







extern int at_quick_exit (void (*__func) (void)) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));







extern int on_exit (void (*__func) (int __status, void *__arg), void *__arg)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));






extern void exit (int __status) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__noreturn__));





extern void quick_exit (int __status) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__noreturn__));







extern void _Exit (int __status) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__noreturn__));






extern char *getenv (const char *__name) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1))) ;

# 577 "/usr/include/stdlib.h" 3 4
extern int putenv (char *__string) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));





extern int setenv (const char *__name, const char *__value, int __replace)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (2)));


extern int unsetenv (const char *__name) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));






extern int clearenv (void) __attribute__ ((__nothrow__ , __leaf__));
# 605 "/usr/include/stdlib.h" 3 4
extern char *mktemp (char *__template) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));
# 619 "/usr/include/stdlib.h" 3 4
extern int mkstemp (char *__template) __attribute__ ((__nonnull__ (1))) ;
# 641 "/usr/include/stdlib.h" 3 4
extern int mkstemps (char *__template, int __suffixlen) __attribute__ ((__nonnull__ (1))) ;
# 662 "/usr/include/stdlib.h" 3 4
extern char *mkdtemp (char *__template) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1))) ;
# 711 "/usr/include/stdlib.h" 3 4





extern int system (const char *__command) ;

# 733 "/usr/include/stdlib.h" 3 4
extern char *realpath (const char *__restrict __name,
         char *__restrict __resolved) __attribute__ ((__nothrow__ , __leaf__)) ;






typedef int (*__compar_fn_t) (const void *, const void *);
# 751 "/usr/include/stdlib.h" 3 4



extern void *bsearch (const void *__key, const void *__base,
        size_t __nmemb, size_t __size, __compar_fn_t __compar)
     __attribute__ ((__nonnull__ (1, 2, 5))) ;



extern void qsort (void *__base, size_t __nmemb, size_t __size,
     __compar_fn_t __compar) __attribute__ ((__nonnull__ (1, 4)));
# 770 "/usr/include/stdlib.h" 3 4
extern int abs (int __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)) ;
extern long int labs (long int __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)) ;



__extension__ extern long long int llabs (long long int __x)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)) ;







extern div_t div (int __numer, int __denom)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)) ;
extern ldiv_t ldiv (long int __numer, long int __denom)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)) ;




__extension__ extern lldiv_t lldiv (long long int __numer,
        long long int __denom)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)) ;

# 807 "/usr/include/stdlib.h" 3 4
extern char *ecvt (double __value, int __ndigit, int *__restrict __decpt,
     int *__restrict __sign) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (3, 4))) ;




extern char *fcvt (double __value, int __ndigit, int *__restrict __decpt,
     int *__restrict __sign) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (3, 4))) ;




extern char *gcvt (double __value, int __ndigit, char *__buf)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (3))) ;




extern char *qecvt (long double __value, int __ndigit,
      int *__restrict __decpt, int *__restrict __sign)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (3, 4))) ;
extern char *qfcvt (long double __value, int __ndigit,
      int *__restrict __decpt, int *__restrict __sign)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (3, 4))) ;
extern char *qgcvt (long double __value, int __ndigit, char *__buf)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (3))) ;




extern int ecvt_r (double __value, int __ndigit, int *__restrict __decpt,
     int *__restrict __sign, char *__restrict __buf,
     size_t __len) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (3, 4, 5)));
extern int fcvt_r (double __value, int __ndigit, int *__restrict __decpt,
     int *__restrict __sign, char *__restrict __buf,
     size_t __len) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (3, 4, 5)));

extern int qecvt_r (long double __value, int __ndigit,
      int *__restrict __decpt, int *__restrict __sign,
      char *__restrict __buf, size_t __len)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (3, 4, 5)));
extern int qfcvt_r (long double __value, int __ndigit,
      int *__restrict __decpt, int *__restrict __sign,
      char *__restrict __buf, size_t __len)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (3, 4, 5)));







extern int mblen (const char *__s, size_t __n) __attribute__ ((__nothrow__ , __leaf__)) ;


extern int mbtowc (wchar_t *__restrict __pwc,
     const char *__restrict __s, size_t __n) __attribute__ ((__nothrow__ , __leaf__)) ;


extern int wctomb (char *__s, wchar_t __wchar) __attribute__ ((__nothrow__ , __leaf__)) ;



extern size_t mbstowcs (wchar_t *__restrict __pwcs,
   const char *__restrict __s, size_t __n) __attribute__ ((__nothrow__ , __leaf__));

extern size_t wcstombs (char *__restrict __s,
   const wchar_t *__restrict __pwcs, size_t __n)
     __attribute__ ((__nothrow__ , __leaf__));








extern int rpmatch (const char *__response) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1))) ;
# 895 "/usr/include/stdlib.h" 3 4
extern int getsubopt (char **__restrict __optionp,
        char *const *__restrict __tokens,
        char **__restrict __valuep)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2, 3))) ;
# 947 "/usr/include/stdlib.h" 3 4
extern int getloadavg (double __loadavg[], int __nelem)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));


# 1 "/usr/include/bits/stdlib-float.h" 1 3 4
# 952 "/usr/include/stdlib.h" 2 3 4
# 964 "/usr/include/stdlib.h" 3 4

# 1873 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h" 2



# 1 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/xmmintrin.h" 1 3 4
# 31 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/xmmintrin.h" 3 4
# 1 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/mmintrin.h" 1 3 4
# 42 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/mmintrin.h" 3 4
typedef int __m64 __attribute__ ((__vector_size__ (8), __may_alias__));


typedef int __m64_u __attribute__ ((__vector_size__ (8), __may_alias__, __aligned__ (1)));


typedef int __v2si __attribute__ ((__vector_size__ (8)));
typedef short __v4hi __attribute__ ((__vector_size__ (8)));
typedef char __v8qi __attribute__ ((__vector_size__ (8)));
typedef long long __v1di __attribute__ ((__vector_size__ (8)));
typedef float __v2sf __attribute__ ((__vector_size__ (8)));


extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_empty (void)
{
  __builtin_ia32_emms ();
}

extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_empty (void)
{
  _mm_empty ();
}


extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtsi32_si64 (int __i)
{
  return (__m64) __builtin_ia32_vec_init_v2si (__i, 0);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_from_int (int __i)
{
  return _mm_cvtsi32_si64 (__i);
}





extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_from_int64 (long long __i)
{
  return (__m64) __i;
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtsi64_m64 (long long __i)
{
  return (__m64) __i;
}


extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtsi64x_si64 (long long __i)
{
  return (__m64) __i;
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_set_pi64x (long long __i)
{
  return (__m64) __i;
}



extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtsi64_si32 (__m64 __i)
{
  return __builtin_ia32_vec_ext_v2si ((__v2si)__i, 0);
}

extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_to_int (__m64 __i)
{
  return _mm_cvtsi64_si32 (__i);
}





extern __inline long long __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_to_int64 (__m64 __i)
{
  return (long long)__i;
}

extern __inline long long __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtm64_si64 (__m64 __i)
{
  return (long long)__i;
}


extern __inline long long __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtsi64_si64x (__m64 __i)
{
  return (long long)__i;
}





extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_packs_pi16 (__m64 __m1, __m64 __m2)
{
  return (__m64) __builtin_ia32_packsswb ((__v4hi)__m1, (__v4hi)__m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_packsswb (__m64 __m1, __m64 __m2)
{
  return _mm_packs_pi16 (__m1, __m2);
}




extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_packs_pi32 (__m64 __m1, __m64 __m2)
{
  return (__m64) __builtin_ia32_packssdw ((__v2si)__m1, (__v2si)__m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_packssdw (__m64 __m1, __m64 __m2)
{
  return _mm_packs_pi32 (__m1, __m2);
}




extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_packs_pu16 (__m64 __m1, __m64 __m2)
{
  return (__m64) __builtin_ia32_packuswb ((__v4hi)__m1, (__v4hi)__m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_packuswb (__m64 __m1, __m64 __m2)
{
  return _mm_packs_pu16 (__m1, __m2);
}



extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_unpackhi_pi8 (__m64 __m1, __m64 __m2)
{
  return (__m64) __builtin_ia32_punpckhbw ((__v8qi)__m1, (__v8qi)__m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_punpckhbw (__m64 __m1, __m64 __m2)
{
  return _mm_unpackhi_pi8 (__m1, __m2);
}



extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_unpackhi_pi16 (__m64 __m1, __m64 __m2)
{
  return (__m64) __builtin_ia32_punpckhwd ((__v4hi)__m1, (__v4hi)__m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_punpckhwd (__m64 __m1, __m64 __m2)
{
  return _mm_unpackhi_pi16 (__m1, __m2);
}



extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_unpackhi_pi32 (__m64 __m1, __m64 __m2)
{
  return (__m64) __builtin_ia32_punpckhdq ((__v2si)__m1, (__v2si)__m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_punpckhdq (__m64 __m1, __m64 __m2)
{
  return _mm_unpackhi_pi32 (__m1, __m2);
}



extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_unpacklo_pi8 (__m64 __m1, __m64 __m2)
{
  return (__m64) __builtin_ia32_punpcklbw ((__v8qi)__m1, (__v8qi)__m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_punpcklbw (__m64 __m1, __m64 __m2)
{
  return _mm_unpacklo_pi8 (__m1, __m2);
}



extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_unpacklo_pi16 (__m64 __m1, __m64 __m2)
{
  return (__m64) __builtin_ia32_punpcklwd ((__v4hi)__m1, (__v4hi)__m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_punpcklwd (__m64 __m1, __m64 __m2)
{
  return _mm_unpacklo_pi16 (__m1, __m2);
}



extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_unpacklo_pi32 (__m64 __m1, __m64 __m2)
{
  return (__m64) __builtin_ia32_punpckldq ((__v2si)__m1, (__v2si)__m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_punpckldq (__m64 __m1, __m64 __m2)
{
  return _mm_unpacklo_pi32 (__m1, __m2);
}


extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_add_pi8 (__m64 __m1, __m64 __m2)
{
  return (__m64) __builtin_ia32_paddb ((__v8qi)__m1, (__v8qi)__m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_paddb (__m64 __m1, __m64 __m2)
{
  return _mm_add_pi8 (__m1, __m2);
}


extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_add_pi16 (__m64 __m1, __m64 __m2)
{
  return (__m64) __builtin_ia32_paddw ((__v4hi)__m1, (__v4hi)__m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_paddw (__m64 __m1, __m64 __m2)
{
  return _mm_add_pi16 (__m1, __m2);
}


extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_add_pi32 (__m64 __m1, __m64 __m2)
{
  return (__m64) __builtin_ia32_paddd ((__v2si)__m1, (__v2si)__m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_paddd (__m64 __m1, __m64 __m2)
{
  return _mm_add_pi32 (__m1, __m2);
}
# 322 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/mmintrin.h" 3 4
extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_add_si64 (__m64 __m1, __m64 __m2)
{
  return (__m64) __builtin_ia32_paddq ((__v1di)__m1, (__v1di)__m2);
}







extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_adds_pi8 (__m64 __m1, __m64 __m2)
{
  return (__m64) __builtin_ia32_paddsb ((__v8qi)__m1, (__v8qi)__m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_paddsb (__m64 __m1, __m64 __m2)
{
  return _mm_adds_pi8 (__m1, __m2);
}



extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_adds_pi16 (__m64 __m1, __m64 __m2)
{
  return (__m64) __builtin_ia32_paddsw ((__v4hi)__m1, (__v4hi)__m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_paddsw (__m64 __m1, __m64 __m2)
{
  return _mm_adds_pi16 (__m1, __m2);
}



extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_adds_pu8 (__m64 __m1, __m64 __m2)
{
  return (__m64) __builtin_ia32_paddusb ((__v8qi)__m1, (__v8qi)__m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_paddusb (__m64 __m1, __m64 __m2)
{
  return _mm_adds_pu8 (__m1, __m2);
}



extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_adds_pu16 (__m64 __m1, __m64 __m2)
{
  return (__m64) __builtin_ia32_paddusw ((__v4hi)__m1, (__v4hi)__m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_paddusw (__m64 __m1, __m64 __m2)
{
  return _mm_adds_pu16 (__m1, __m2);
}


extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_sub_pi8 (__m64 __m1, __m64 __m2)
{
  return (__m64) __builtin_ia32_psubb ((__v8qi)__m1, (__v8qi)__m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_psubb (__m64 __m1, __m64 __m2)
{
  return _mm_sub_pi8 (__m1, __m2);
}


extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_sub_pi16 (__m64 __m1, __m64 __m2)
{
  return (__m64) __builtin_ia32_psubw ((__v4hi)__m1, (__v4hi)__m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_psubw (__m64 __m1, __m64 __m2)
{
  return _mm_sub_pi16 (__m1, __m2);
}


extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_sub_pi32 (__m64 __m1, __m64 __m2)
{
  return (__m64) __builtin_ia32_psubd ((__v2si)__m1, (__v2si)__m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_psubd (__m64 __m1, __m64 __m2)
{
  return _mm_sub_pi32 (__m1, __m2);
}
# 434 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/mmintrin.h" 3 4
extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_sub_si64 (__m64 __m1, __m64 __m2)
{
  return (__m64) __builtin_ia32_psubq ((__v1di)__m1, (__v1di)__m2);
}







extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_subs_pi8 (__m64 __m1, __m64 __m2)
{
  return (__m64) __builtin_ia32_psubsb ((__v8qi)__m1, (__v8qi)__m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_psubsb (__m64 __m1, __m64 __m2)
{
  return _mm_subs_pi8 (__m1, __m2);
}



extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_subs_pi16 (__m64 __m1, __m64 __m2)
{
  return (__m64) __builtin_ia32_psubsw ((__v4hi)__m1, (__v4hi)__m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_psubsw (__m64 __m1, __m64 __m2)
{
  return _mm_subs_pi16 (__m1, __m2);
}



extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_subs_pu8 (__m64 __m1, __m64 __m2)
{
  return (__m64) __builtin_ia32_psubusb ((__v8qi)__m1, (__v8qi)__m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_psubusb (__m64 __m1, __m64 __m2)
{
  return _mm_subs_pu8 (__m1, __m2);
}



extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_subs_pu16 (__m64 __m1, __m64 __m2)
{
  return (__m64) __builtin_ia32_psubusw ((__v4hi)__m1, (__v4hi)__m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_psubusw (__m64 __m1, __m64 __m2)
{
  return _mm_subs_pu16 (__m1, __m2);
}




extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_madd_pi16 (__m64 __m1, __m64 __m2)
{
  return (__m64) __builtin_ia32_pmaddwd ((__v4hi)__m1, (__v4hi)__m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_pmaddwd (__m64 __m1, __m64 __m2)
{
  return _mm_madd_pi16 (__m1, __m2);
}



extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_mulhi_pi16 (__m64 __m1, __m64 __m2)
{
  return (__m64) __builtin_ia32_pmulhw ((__v4hi)__m1, (__v4hi)__m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_pmulhw (__m64 __m1, __m64 __m2)
{
  return _mm_mulhi_pi16 (__m1, __m2);
}



extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_mullo_pi16 (__m64 __m1, __m64 __m2)
{
  return (__m64) __builtin_ia32_pmullw ((__v4hi)__m1, (__v4hi)__m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_pmullw (__m64 __m1, __m64 __m2)
{
  return _mm_mullo_pi16 (__m1, __m2);
}


extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_sll_pi16 (__m64 __m, __m64 __count)
{
  return (__m64) __builtin_ia32_psllw ((__v4hi)__m, (__v4hi)__count);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_psllw (__m64 __m, __m64 __count)
{
  return _mm_sll_pi16 (__m, __count);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_slli_pi16 (__m64 __m, int __count)
{
  return (__m64) __builtin_ia32_psllwi ((__v4hi)__m, __count);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_psllwi (__m64 __m, int __count)
{
  return _mm_slli_pi16 (__m, __count);
}


extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_sll_pi32 (__m64 __m, __m64 __count)
{
  return (__m64) __builtin_ia32_pslld ((__v2si)__m, (__v2si)__count);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_pslld (__m64 __m, __m64 __count)
{
  return _mm_sll_pi32 (__m, __count);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_slli_pi32 (__m64 __m, int __count)
{
  return (__m64) __builtin_ia32_pslldi ((__v2si)__m, __count);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_pslldi (__m64 __m, int __count)
{
  return _mm_slli_pi32 (__m, __count);
}


extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_sll_si64 (__m64 __m, __m64 __count)
{
  return (__m64) __builtin_ia32_psllq ((__v1di)__m, (__v1di)__count);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_psllq (__m64 __m, __m64 __count)
{
  return _mm_sll_si64 (__m, __count);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_slli_si64 (__m64 __m, int __count)
{
  return (__m64) __builtin_ia32_psllqi ((__v1di)__m, __count);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_psllqi (__m64 __m, int __count)
{
  return _mm_slli_si64 (__m, __count);
}


extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_sra_pi16 (__m64 __m, __m64 __count)
{
  return (__m64) __builtin_ia32_psraw ((__v4hi)__m, (__v4hi)__count);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_psraw (__m64 __m, __m64 __count)
{
  return _mm_sra_pi16 (__m, __count);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_srai_pi16 (__m64 __m, int __count)
{
  return (__m64) __builtin_ia32_psrawi ((__v4hi)__m, __count);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_psrawi (__m64 __m, int __count)
{
  return _mm_srai_pi16 (__m, __count);
}


extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_sra_pi32 (__m64 __m, __m64 __count)
{
  return (__m64) __builtin_ia32_psrad ((__v2si)__m, (__v2si)__count);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_psrad (__m64 __m, __m64 __count)
{
  return _mm_sra_pi32 (__m, __count);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_srai_pi32 (__m64 __m, int __count)
{
  return (__m64) __builtin_ia32_psradi ((__v2si)__m, __count);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_psradi (__m64 __m, int __count)
{
  return _mm_srai_pi32 (__m, __count);
}


extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_srl_pi16 (__m64 __m, __m64 __count)
{
  return (__m64) __builtin_ia32_psrlw ((__v4hi)__m, (__v4hi)__count);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_psrlw (__m64 __m, __m64 __count)
{
  return _mm_srl_pi16 (__m, __count);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_srli_pi16 (__m64 __m, int __count)
{
  return (__m64) __builtin_ia32_psrlwi ((__v4hi)__m, __count);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_psrlwi (__m64 __m, int __count)
{
  return _mm_srli_pi16 (__m, __count);
}


extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_srl_pi32 (__m64 __m, __m64 __count)
{
  return (__m64) __builtin_ia32_psrld ((__v2si)__m, (__v2si)__count);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_psrld (__m64 __m, __m64 __count)
{
  return _mm_srl_pi32 (__m, __count);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_srli_pi32 (__m64 __m, int __count)
{
  return (__m64) __builtin_ia32_psrldi ((__v2si)__m, __count);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_psrldi (__m64 __m, int __count)
{
  return _mm_srli_pi32 (__m, __count);
}


extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_srl_si64 (__m64 __m, __m64 __count)
{
  return (__m64) __builtin_ia32_psrlq ((__v1di)__m, (__v1di)__count);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_psrlq (__m64 __m, __m64 __count)
{
  return _mm_srl_si64 (__m, __count);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_srli_si64 (__m64 __m, int __count)
{
  return (__m64) __builtin_ia32_psrlqi ((__v1di)__m, __count);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_psrlqi (__m64 __m, int __count)
{
  return _mm_srli_si64 (__m, __count);
}


extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_and_si64 (__m64 __m1, __m64 __m2)
{
  return __builtin_ia32_pand (__m1, __m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_pand (__m64 __m1, __m64 __m2)
{
  return _mm_and_si64 (__m1, __m2);
}



extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_andnot_si64 (__m64 __m1, __m64 __m2)
{
  return __builtin_ia32_pandn (__m1, __m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_pandn (__m64 __m1, __m64 __m2)
{
  return _mm_andnot_si64 (__m1, __m2);
}


extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_or_si64 (__m64 __m1, __m64 __m2)
{
  return __builtin_ia32_por (__m1, __m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_por (__m64 __m1, __m64 __m2)
{
  return _mm_or_si64 (__m1, __m2);
}


extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_xor_si64 (__m64 __m1, __m64 __m2)
{
  return __builtin_ia32_pxor (__m1, __m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_pxor (__m64 __m1, __m64 __m2)
{
  return _mm_xor_si64 (__m1, __m2);
}



extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpeq_pi8 (__m64 __m1, __m64 __m2)
{
  return (__m64) __builtin_ia32_pcmpeqb ((__v8qi)__m1, (__v8qi)__m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_pcmpeqb (__m64 __m1, __m64 __m2)
{
  return _mm_cmpeq_pi8 (__m1, __m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpgt_pi8 (__m64 __m1, __m64 __m2)
{
  return (__m64) __builtin_ia32_pcmpgtb ((__v8qi)__m1, (__v8qi)__m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_pcmpgtb (__m64 __m1, __m64 __m2)
{
  return _mm_cmpgt_pi8 (__m1, __m2);
}



extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpeq_pi16 (__m64 __m1, __m64 __m2)
{
  return (__m64) __builtin_ia32_pcmpeqw ((__v4hi)__m1, (__v4hi)__m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_pcmpeqw (__m64 __m1, __m64 __m2)
{
  return _mm_cmpeq_pi16 (__m1, __m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpgt_pi16 (__m64 __m1, __m64 __m2)
{
  return (__m64) __builtin_ia32_pcmpgtw ((__v4hi)__m1, (__v4hi)__m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_pcmpgtw (__m64 __m1, __m64 __m2)
{
  return _mm_cmpgt_pi16 (__m1, __m2);
}



extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpeq_pi32 (__m64 __m1, __m64 __m2)
{
  return (__m64) __builtin_ia32_pcmpeqd ((__v2si)__m1, (__v2si)__m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_pcmpeqd (__m64 __m1, __m64 __m2)
{
  return _mm_cmpeq_pi32 (__m1, __m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpgt_pi32 (__m64 __m1, __m64 __m2)
{
  return (__m64) __builtin_ia32_pcmpgtd ((__v2si)__m1, (__v2si)__m2);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_pcmpgtd (__m64 __m1, __m64 __m2)
{
  return _mm_cmpgt_pi32 (__m1, __m2);
}


extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_setzero_si64 (void)
{
  return (__m64)0LL;
}


extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_set_pi32 (int __i1, int __i0)
{
  return (__m64) __builtin_ia32_vec_init_v2si (__i0, __i1);
}


extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_set_pi16 (short __w3, short __w2, short __w1, short __w0)
{
  return (__m64) __builtin_ia32_vec_init_v4hi (__w0, __w1, __w2, __w3);
}


extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_set_pi8 (char __b7, char __b6, char __b5, char __b4,
      char __b3, char __b2, char __b1, char __b0)
{
  return (__m64) __builtin_ia32_vec_init_v8qi (__b0, __b1, __b2, __b3,
            __b4, __b5, __b6, __b7);
}


extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_setr_pi32 (int __i0, int __i1)
{
  return _mm_set_pi32 (__i1, __i0);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_setr_pi16 (short __w0, short __w1, short __w2, short __w3)
{
  return _mm_set_pi16 (__w3, __w2, __w1, __w0);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_setr_pi8 (char __b0, char __b1, char __b2, char __b3,
       char __b4, char __b5, char __b6, char __b7)
{
  return _mm_set_pi8 (__b7, __b6, __b5, __b4, __b3, __b2, __b1, __b0);
}


extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_set1_pi32 (int __i)
{
  return _mm_set_pi32 (__i, __i);
}


extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_set1_pi16 (short __w)
{
  return _mm_set_pi16 (__w, __w, __w, __w);
}


extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_set1_pi8 (char __b)
{
  return _mm_set_pi8 (__b, __b, __b, __b, __b, __b, __b, __b);
}
# 32 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/xmmintrin.h" 2 3 4


# 1 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/mm_malloc.h" 1 3 4
# 32 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/mm_malloc.h" 3 4
extern int posix_memalign (void **, size_t, size_t);




static __inline void *
_mm_malloc (size_t __size, size_t __alignment)
{
  void *__ptr;
  if (__alignment == 1)
    return malloc (__size);
  if (__alignment == 2 || (sizeof (void *) == 8 && __alignment == 4))
    __alignment = sizeof (void *);
  if (posix_memalign (&__ptr, __alignment, __size) == 0)
    return __ptr;
  else
    return ((void *)0);
}

static __inline void
_mm_free (void *__ptr)
{
  free (__ptr);
}
# 35 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/xmmintrin.h" 2 3 4


enum _mm_hint
{

  _MM_HINT_ET0 = 7,
  _MM_HINT_ET1 = 6,
  _MM_HINT_T0 = 3,
  _MM_HINT_T1 = 2,
  _MM_HINT_T2 = 1,
  _MM_HINT_NTA = 0
};
# 69 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/xmmintrin.h" 3 4
typedef float __m128 __attribute__ ((__vector_size__ (16), __may_alias__));


typedef float __m128_u __attribute__ ((__vector_size__ (16), __may_alias__, __aligned__ (1)));


typedef float __v4sf __attribute__ ((__vector_size__ (16)));
# 109 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/xmmintrin.h" 3 4
extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_undefined_ps (void)
{
  __m128 __Y = __Y;
  return __Y;
}


extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_setzero_ps (void)
{
  return __extension__ (__m128){ 0.0f, 0.0f, 0.0f, 0.0f };
}





extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_add_ss (__m128 __A, __m128 __B)
{
  return (__m128) __builtin_ia32_addss ((__v4sf)__A, (__v4sf)__B);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_sub_ss (__m128 __A, __m128 __B)
{
  return (__m128) __builtin_ia32_subss ((__v4sf)__A, (__v4sf)__B);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_mul_ss (__m128 __A, __m128 __B)
{
  return (__m128) __builtin_ia32_mulss ((__v4sf)__A, (__v4sf)__B);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_div_ss (__m128 __A, __m128 __B)
{
  return (__m128) __builtin_ia32_divss ((__v4sf)__A, (__v4sf)__B);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_sqrt_ss (__m128 __A)
{
  return (__m128) __builtin_ia32_sqrtss ((__v4sf)__A);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_rcp_ss (__m128 __A)
{
  return (__m128) __builtin_ia32_rcpss ((__v4sf)__A);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_rsqrt_ss (__m128 __A)
{
  return (__m128) __builtin_ia32_rsqrtss ((__v4sf)__A);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_min_ss (__m128 __A, __m128 __B)
{
  return (__m128) __builtin_ia32_minss ((__v4sf)__A, (__v4sf)__B);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_max_ss (__m128 __A, __m128 __B)
{
  return (__m128) __builtin_ia32_maxss ((__v4sf)__A, (__v4sf)__B);
}



extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_add_ps (__m128 __A, __m128 __B)
{
  return (__m128) ((__v4sf)__A + (__v4sf)__B);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_sub_ps (__m128 __A, __m128 __B)
{
  return (__m128) ((__v4sf)__A - (__v4sf)__B);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_mul_ps (__m128 __A, __m128 __B)
{
  return (__m128) ((__v4sf)__A * (__v4sf)__B);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_div_ps (__m128 __A, __m128 __B)
{
  return (__m128) ((__v4sf)__A / (__v4sf)__B);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_sqrt_ps (__m128 __A)
{
  return (__m128) __builtin_ia32_sqrtps ((__v4sf)__A);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_rcp_ps (__m128 __A)
{
  return (__m128) __builtin_ia32_rcpps ((__v4sf)__A);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_rsqrt_ps (__m128 __A)
{
  return (__m128) __builtin_ia32_rsqrtps ((__v4sf)__A);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_min_ps (__m128 __A, __m128 __B)
{
  return (__m128) __builtin_ia32_minps ((__v4sf)__A, (__v4sf)__B);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_max_ps (__m128 __A, __m128 __B)
{
  return (__m128) __builtin_ia32_maxps ((__v4sf)__A, (__v4sf)__B);
}



extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_and_ps (__m128 __A, __m128 __B)
{
  return __builtin_ia32_andps (__A, __B);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_andnot_ps (__m128 __A, __m128 __B)
{
  return __builtin_ia32_andnps (__A, __B);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_or_ps (__m128 __A, __m128 __B)
{
  return __builtin_ia32_orps (__A, __B);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_xor_ps (__m128 __A, __m128 __B)
{
  return __builtin_ia32_xorps (__A, __B);
}





extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpeq_ss (__m128 __A, __m128 __B)
{
  return (__m128) __builtin_ia32_cmpeqss ((__v4sf)__A, (__v4sf)__B);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmplt_ss (__m128 __A, __m128 __B)
{
  return (__m128) __builtin_ia32_cmpltss ((__v4sf)__A, (__v4sf)__B);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmple_ss (__m128 __A, __m128 __B)
{
  return (__m128) __builtin_ia32_cmpless ((__v4sf)__A, (__v4sf)__B);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpgt_ss (__m128 __A, __m128 __B)
{
  return (__m128) __builtin_ia32_movss ((__v4sf) __A,
     (__v4sf)
     __builtin_ia32_cmpltss ((__v4sf) __B,
        (__v4sf)
        __A));
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpge_ss (__m128 __A, __m128 __B)
{
  return (__m128) __builtin_ia32_movss ((__v4sf) __A,
     (__v4sf)
     __builtin_ia32_cmpless ((__v4sf) __B,
        (__v4sf)
        __A));
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpneq_ss (__m128 __A, __m128 __B)
{
  return (__m128) __builtin_ia32_cmpneqss ((__v4sf)__A, (__v4sf)__B);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpnlt_ss (__m128 __A, __m128 __B)
{
  return (__m128) __builtin_ia32_cmpnltss ((__v4sf)__A, (__v4sf)__B);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpnle_ss (__m128 __A, __m128 __B)
{
  return (__m128) __builtin_ia32_cmpnless ((__v4sf)__A, (__v4sf)__B);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpngt_ss (__m128 __A, __m128 __B)
{
  return (__m128) __builtin_ia32_movss ((__v4sf) __A,
     (__v4sf)
     __builtin_ia32_cmpnltss ((__v4sf) __B,
         (__v4sf)
         __A));
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpnge_ss (__m128 __A, __m128 __B)
{
  return (__m128) __builtin_ia32_movss ((__v4sf) __A,
     (__v4sf)
     __builtin_ia32_cmpnless ((__v4sf) __B,
         (__v4sf)
         __A));
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpord_ss (__m128 __A, __m128 __B)
{
  return (__m128) __builtin_ia32_cmpordss ((__v4sf)__A, (__v4sf)__B);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpunord_ss (__m128 __A, __m128 __B)
{
  return (__m128) __builtin_ia32_cmpunordss ((__v4sf)__A, (__v4sf)__B);
}





extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpeq_ps (__m128 __A, __m128 __B)
{
  return (__m128) __builtin_ia32_cmpeqps ((__v4sf)__A, (__v4sf)__B);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmplt_ps (__m128 __A, __m128 __B)
{
  return (__m128) __builtin_ia32_cmpltps ((__v4sf)__A, (__v4sf)__B);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmple_ps (__m128 __A, __m128 __B)
{
  return (__m128) __builtin_ia32_cmpleps ((__v4sf)__A, (__v4sf)__B);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpgt_ps (__m128 __A, __m128 __B)
{
  return (__m128) __builtin_ia32_cmpgtps ((__v4sf)__A, (__v4sf)__B);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpge_ps (__m128 __A, __m128 __B)
{
  return (__m128) __builtin_ia32_cmpgeps ((__v4sf)__A, (__v4sf)__B);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpneq_ps (__m128 __A, __m128 __B)
{
  return (__m128) __builtin_ia32_cmpneqps ((__v4sf)__A, (__v4sf)__B);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpnlt_ps (__m128 __A, __m128 __B)
{
  return (__m128) __builtin_ia32_cmpnltps ((__v4sf)__A, (__v4sf)__B);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpnle_ps (__m128 __A, __m128 __B)
{
  return (__m128) __builtin_ia32_cmpnleps ((__v4sf)__A, (__v4sf)__B);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpngt_ps (__m128 __A, __m128 __B)
{
  return (__m128) __builtin_ia32_cmpngtps ((__v4sf)__A, (__v4sf)__B);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpnge_ps (__m128 __A, __m128 __B)
{
  return (__m128) __builtin_ia32_cmpngeps ((__v4sf)__A, (__v4sf)__B);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpord_ps (__m128 __A, __m128 __B)
{
  return (__m128) __builtin_ia32_cmpordps ((__v4sf)__A, (__v4sf)__B);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpunord_ps (__m128 __A, __m128 __B)
{
  return (__m128) __builtin_ia32_cmpunordps ((__v4sf)__A, (__v4sf)__B);
}




extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_comieq_ss (__m128 __A, __m128 __B)
{
  return __builtin_ia32_comieq ((__v4sf)__A, (__v4sf)__B);
}

extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_comilt_ss (__m128 __A, __m128 __B)
{
  return __builtin_ia32_comilt ((__v4sf)__A, (__v4sf)__B);
}

extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_comile_ss (__m128 __A, __m128 __B)
{
  return __builtin_ia32_comile ((__v4sf)__A, (__v4sf)__B);
}

extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_comigt_ss (__m128 __A, __m128 __B)
{
  return __builtin_ia32_comigt ((__v4sf)__A, (__v4sf)__B);
}

extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_comige_ss (__m128 __A, __m128 __B)
{
  return __builtin_ia32_comige ((__v4sf)__A, (__v4sf)__B);
}

extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_comineq_ss (__m128 __A, __m128 __B)
{
  return __builtin_ia32_comineq ((__v4sf)__A, (__v4sf)__B);
}

extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_ucomieq_ss (__m128 __A, __m128 __B)
{
  return __builtin_ia32_ucomieq ((__v4sf)__A, (__v4sf)__B);
}

extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_ucomilt_ss (__m128 __A, __m128 __B)
{
  return __builtin_ia32_ucomilt ((__v4sf)__A, (__v4sf)__B);
}

extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_ucomile_ss (__m128 __A, __m128 __B)
{
  return __builtin_ia32_ucomile ((__v4sf)__A, (__v4sf)__B);
}

extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_ucomigt_ss (__m128 __A, __m128 __B)
{
  return __builtin_ia32_ucomigt ((__v4sf)__A, (__v4sf)__B);
}

extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_ucomige_ss (__m128 __A, __m128 __B)
{
  return __builtin_ia32_ucomige ((__v4sf)__A, (__v4sf)__B);
}

extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_ucomineq_ss (__m128 __A, __m128 __B)
{
  return __builtin_ia32_ucomineq ((__v4sf)__A, (__v4sf)__B);
}



extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtss_si32 (__m128 __A)
{
  return __builtin_ia32_cvtss2si ((__v4sf) __A);
}

extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvt_ss2si (__m128 __A)
{
  return _mm_cvtss_si32 (__A);
}






extern __inline long long __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtss_si64 (__m128 __A)
{
  return __builtin_ia32_cvtss2si64 ((__v4sf) __A);
}


extern __inline long long __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtss_si64x (__m128 __A)
{
  return __builtin_ia32_cvtss2si64 ((__v4sf) __A);
}




extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtps_pi32 (__m128 __A)
{
  return (__m64) __builtin_ia32_cvtps2pi ((__v4sf) __A);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvt_ps2pi (__m128 __A)
{
  return _mm_cvtps_pi32 (__A);
}


extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvttss_si32 (__m128 __A)
{
  return __builtin_ia32_cvttss2si ((__v4sf) __A);
}

extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtt_ss2si (__m128 __A)
{
  return _mm_cvttss_si32 (__A);
}





extern __inline long long __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvttss_si64 (__m128 __A)
{
  return __builtin_ia32_cvttss2si64 ((__v4sf) __A);
}


extern __inline long long __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvttss_si64x (__m128 __A)
{
  return __builtin_ia32_cvttss2si64 ((__v4sf) __A);
}




extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvttps_pi32 (__m128 __A)
{
  return (__m64) __builtin_ia32_cvttps2pi ((__v4sf) __A);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtt_ps2pi (__m128 __A)
{
  return _mm_cvttps_pi32 (__A);
}


extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtsi32_ss (__m128 __A, int __B)
{
  return (__m128) __builtin_ia32_cvtsi2ss ((__v4sf) __A, __B);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvt_si2ss (__m128 __A, int __B)
{
  return _mm_cvtsi32_ss (__A, __B);
}





extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtsi64_ss (__m128 __A, long long __B)
{
  return (__m128) __builtin_ia32_cvtsi642ss ((__v4sf) __A, __B);
}


extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtsi64x_ss (__m128 __A, long long __B)
{
  return (__m128) __builtin_ia32_cvtsi642ss ((__v4sf) __A, __B);
}




extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtpi32_ps (__m128 __A, __m64 __B)
{
  return (__m128) __builtin_ia32_cvtpi2ps ((__v4sf) __A, (__v2si)__B);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvt_pi2ps (__m128 __A, __m64 __B)
{
  return _mm_cvtpi32_ps (__A, __B);
}


extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtpi16_ps (__m64 __A)
{
  __v4hi __sign;
  __v2si __hisi, __losi;
  __v4sf __zero, __ra, __rb;




  __sign = __builtin_ia32_pcmpgtw ((__v4hi)0LL, (__v4hi)__A);


  __losi = (__v2si) __builtin_ia32_punpcklwd ((__v4hi)__A, __sign);
  __hisi = (__v2si) __builtin_ia32_punpckhwd ((__v4hi)__A, __sign);


  __zero = (__v4sf) _mm_setzero_ps ();
  __ra = __builtin_ia32_cvtpi2ps (__zero, __losi);
  __rb = __builtin_ia32_cvtpi2ps (__ra, __hisi);

  return (__m128) __builtin_ia32_movlhps (__ra, __rb);
}


extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtpu16_ps (__m64 __A)
{
  __v2si __hisi, __losi;
  __v4sf __zero, __ra, __rb;


  __losi = (__v2si) __builtin_ia32_punpcklwd ((__v4hi)__A, (__v4hi)0LL);
  __hisi = (__v2si) __builtin_ia32_punpckhwd ((__v4hi)__A, (__v4hi)0LL);


  __zero = (__v4sf) _mm_setzero_ps ();
  __ra = __builtin_ia32_cvtpi2ps (__zero, __losi);
  __rb = __builtin_ia32_cvtpi2ps (__ra, __hisi);

  return (__m128) __builtin_ia32_movlhps (__ra, __rb);
}


extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtpi8_ps (__m64 __A)
{
  __v8qi __sign;




  __sign = __builtin_ia32_pcmpgtb ((__v8qi)0LL, (__v8qi)__A);


  __A = (__m64) __builtin_ia32_punpcklbw ((__v8qi)__A, __sign);

  return _mm_cvtpi16_ps(__A);
}


extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtpu8_ps(__m64 __A)
{
  __A = (__m64) __builtin_ia32_punpcklbw ((__v8qi)__A, (__v8qi)0LL);
  return _mm_cvtpu16_ps(__A);
}


extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtpi32x2_ps(__m64 __A, __m64 __B)
{
  __v4sf __zero = (__v4sf) _mm_setzero_ps ();
  __v4sf __sfa = __builtin_ia32_cvtpi2ps (__zero, (__v2si)__A);
  __v4sf __sfb = __builtin_ia32_cvtpi2ps (__sfa, (__v2si)__B);
  return (__m128) __builtin_ia32_movlhps (__sfa, __sfb);
}


extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtps_pi16(__m128 __A)
{
  __v4sf __hisf = (__v4sf)__A;
  __v4sf __losf = __builtin_ia32_movhlps (__hisf, __hisf);
  __v2si __hisi = __builtin_ia32_cvtps2pi (__hisf);
  __v2si __losi = __builtin_ia32_cvtps2pi (__losf);
  return (__m64) __builtin_ia32_packssdw (__hisi, __losi);
}


extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtps_pi8(__m128 __A)
{
  __v4hi __tmp = (__v4hi) _mm_cvtps_pi16 (__A);
  return (__m64) __builtin_ia32_packsswb (__tmp, (__v4hi)0LL);
}
# 755 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/xmmintrin.h" 3 4
extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_unpackhi_ps (__m128 __A, __m128 __B)
{
  return (__m128) __builtin_ia32_unpckhps ((__v4sf)__A, (__v4sf)__B);
}


extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_unpacklo_ps (__m128 __A, __m128 __B)
{
  return (__m128) __builtin_ia32_unpcklps ((__v4sf)__A, (__v4sf)__B);
}



extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_loadh_pi (__m128 __A, __m64 const *__P)
{
  return (__m128) __builtin_ia32_loadhps ((__v4sf)__A, (const __v2sf *)__P);
}


extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_storeh_pi (__m64 *__P, __m128 __A)
{
  __builtin_ia32_storehps ((__v2sf *)__P, (__v4sf)__A);
}


extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_movehl_ps (__m128 __A, __m128 __B)
{
  return (__m128) __builtin_ia32_movhlps ((__v4sf)__A, (__v4sf)__B);
}


extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_movelh_ps (__m128 __A, __m128 __B)
{
  return (__m128) __builtin_ia32_movlhps ((__v4sf)__A, (__v4sf)__B);
}



extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_loadl_pi (__m128 __A, __m64 const *__P)
{
  return (__m128) __builtin_ia32_loadlps ((__v4sf)__A, (const __v2sf *)__P);
}


extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_storel_pi (__m64 *__P, __m128 __A)
{
  __builtin_ia32_storelps ((__v2sf *)__P, (__v4sf)__A);
}


extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_movemask_ps (__m128 __A)
{
  return __builtin_ia32_movmskps ((__v4sf)__A);
}


extern __inline unsigned int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_getcsr (void)
{
  return __builtin_ia32_stmxcsr ();
}


extern __inline unsigned int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_MM_GET_EXCEPTION_STATE (void)
{
  return _mm_getcsr() & 0x003f;
}

extern __inline unsigned int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_MM_GET_EXCEPTION_MASK (void)
{
  return _mm_getcsr() & 0x1f80;
}

extern __inline unsigned int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_MM_GET_ROUNDING_MODE (void)
{
  return _mm_getcsr() & 0x6000;
}

extern __inline unsigned int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_MM_GET_FLUSH_ZERO_MODE (void)
{
  return _mm_getcsr() & 0x8000;
}


extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_setcsr (unsigned int __I)
{
  __builtin_ia32_ldmxcsr (__I);
}


extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_MM_SET_EXCEPTION_STATE(unsigned int __mask)
{
  _mm_setcsr((_mm_getcsr() & ~0x003f) | __mask);
}

extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_MM_SET_EXCEPTION_MASK (unsigned int __mask)
{
  _mm_setcsr((_mm_getcsr() & ~0x1f80) | __mask);
}

extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_MM_SET_ROUNDING_MODE (unsigned int __mode)
{
  _mm_setcsr((_mm_getcsr() & ~0x6000) | __mode);
}

extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_MM_SET_FLUSH_ZERO_MODE (unsigned int __mode)
{
  _mm_setcsr((_mm_getcsr() & ~0x8000) | __mode);
}


extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_set_ss (float __F)
{
  return __extension__ (__m128)(__v4sf){ __F, 0.0f, 0.0f, 0.0f };
}


extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_set1_ps (float __F)
{
  return __extension__ (__m128)(__v4sf){ __F, __F, __F, __F };
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_set_ps1 (float __F)
{
  return _mm_set1_ps (__F);
}


extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_load_ss (float const *__P)
{
  return _mm_set_ss (*__P);
}


extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_load1_ps (float const *__P)
{
  return _mm_set1_ps (*__P);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_load_ps1 (float const *__P)
{
  return _mm_load1_ps (__P);
}


extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_load_ps (float const *__P)
{
  return *(__m128 *)__P;
}


extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_loadu_ps (float const *__P)
{
  return *(__m128_u *)__P;
}


extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_loadr_ps (float const *__P)
{
  __v4sf __tmp = *(__v4sf *)__P;
  return (__m128) __builtin_ia32_shufps (__tmp, __tmp, (((0) << 6) | ((1) << 4) | ((2) << 2) | (3)));
}


extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_set_ps (const float __Z, const float __Y, const float __X, const float __W)
{
  return __extension__ (__m128)(__v4sf){ __W, __X, __Y, __Z };
}


extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_setr_ps (float __Z, float __Y, float __X, float __W)
{
  return __extension__ (__m128)(__v4sf){ __Z, __Y, __X, __W };
}


extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_store_ss (float *__P, __m128 __A)
{
  *__P = ((__v4sf)__A)[0];
}

extern __inline float __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtss_f32 (__m128 __A)
{
  return ((__v4sf)__A)[0];
}


extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_store_ps (float *__P, __m128 __A)
{
  *(__m128 *)__P = __A;
}


extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_storeu_ps (float *__P, __m128 __A)
{
  *(__m128_u *)__P = __A;
}


extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_store1_ps (float *__P, __m128 __A)
{
  __v4sf __va = (__v4sf)__A;
  __v4sf __tmp = __builtin_ia32_shufps (__va, __va, (((0) << 6) | ((0) << 4) | ((0) << 2) | (0)));
  _mm_storeu_ps (__P, __tmp);
}

extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_store_ps1 (float *__P, __m128 __A)
{
  _mm_store1_ps (__P, __A);
}


extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_storer_ps (float *__P, __m128 __A)
{
  __v4sf __va = (__v4sf)__A;
  __v4sf __tmp = __builtin_ia32_shufps (__va, __va, (((0) << 6) | ((1) << 4) | ((2) << 2) | (3)));
  _mm_store_ps (__P, __tmp);
}


extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_move_ss (__m128 __A, __m128 __B)
{
  return (__m128) __builtin_ia32_movss ((__v4sf)__A, (__v4sf)__B);
}
# 1060 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/xmmintrin.h" 3 4
extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_max_pi16 (__m64 __A, __m64 __B)
{
  return (__m64) __builtin_ia32_pmaxsw ((__v4hi)__A, (__v4hi)__B);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_pmaxsw (__m64 __A, __m64 __B)
{
  return _mm_max_pi16 (__A, __B);
}


extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_max_pu8 (__m64 __A, __m64 __B)
{
  return (__m64) __builtin_ia32_pmaxub ((__v8qi)__A, (__v8qi)__B);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_pmaxub (__m64 __A, __m64 __B)
{
  return _mm_max_pu8 (__A, __B);
}


extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_min_pi16 (__m64 __A, __m64 __B)
{
  return (__m64) __builtin_ia32_pminsw ((__v4hi)__A, (__v4hi)__B);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_pminsw (__m64 __A, __m64 __B)
{
  return _mm_min_pi16 (__A, __B);
}


extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_min_pu8 (__m64 __A, __m64 __B)
{
  return (__m64) __builtin_ia32_pminub ((__v8qi)__A, (__v8qi)__B);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_pminub (__m64 __A, __m64 __B)
{
  return _mm_min_pu8 (__A, __B);
}


extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_movemask_pi8 (__m64 __A)
{
  return __builtin_ia32_pmovmskb ((__v8qi)__A);
}

extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_pmovmskb (__m64 __A)
{
  return _mm_movemask_pi8 (__A);
}



extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_mulhi_pu16 (__m64 __A, __m64 __B)
{
  return (__m64) __builtin_ia32_pmulhuw ((__v4hi)__A, (__v4hi)__B);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_pmulhuw (__m64 __A, __m64 __B)
{
  return _mm_mulhi_pu16 (__A, __B);
}
# 1162 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/xmmintrin.h" 3 4
extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_maskmove_si64 (__m64 __A, __m64 __N, char *__P)
{
  __builtin_ia32_maskmovq ((__v8qi)__A, (__v8qi)__N, __P);
}

extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_maskmovq (__m64 __A, __m64 __N, char *__P)
{
  _mm_maskmove_si64 (__A, __N, __P);
}


extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_avg_pu8 (__m64 __A, __m64 __B)
{
  return (__m64) __builtin_ia32_pavgb ((__v8qi)__A, (__v8qi)__B);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_pavgb (__m64 __A, __m64 __B)
{
  return _mm_avg_pu8 (__A, __B);
}


extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_avg_pu16 (__m64 __A, __m64 __B)
{
  return (__m64) __builtin_ia32_pavgw ((__v4hi)__A, (__v4hi)__B);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_pavgw (__m64 __A, __m64 __B)
{
  return _mm_avg_pu16 (__A, __B);
}




extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_sad_pu8 (__m64 __A, __m64 __B)
{
  return (__m64) __builtin_ia32_psadbw ((__v8qi)__A, (__v8qi)__B);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_m_psadbw (__m64 __A, __m64 __B)
{
  return _mm_sad_pu8 (__A, __B);
}


extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_stream_pi (__m64 *__P, __m64 __A)
{
  __builtin_ia32_movntq ((unsigned long long *)__P, (unsigned long long)__A);
}


extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_stream_ps (float *__P, __m128 __A)
{
  __builtin_ia32_movntps (__P, (__v4sf)__A);
}



extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_sfence (void)
{
  __builtin_ia32_sfence ();
}
# 1252 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/xmmintrin.h" 3 4
# 1 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/emmintrin.h" 1 3 4
# 31 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/emmintrin.h" 3 4
# 1 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/xmmintrin.h" 1 3 4
# 32 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/emmintrin.h" 2 3 4
# 40 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/emmintrin.h" 3 4
typedef double __v2df __attribute__ ((__vector_size__ (16)));
typedef long long __v2di __attribute__ ((__vector_size__ (16)));
typedef unsigned long long __v2du __attribute__ ((__vector_size__ (16)));
typedef int __v4si __attribute__ ((__vector_size__ (16)));
typedef unsigned int __v4su __attribute__ ((__vector_size__ (16)));
typedef short __v8hi __attribute__ ((__vector_size__ (16)));
typedef unsigned short __v8hu __attribute__ ((__vector_size__ (16)));
typedef char __v16qi __attribute__ ((__vector_size__ (16)));
typedef signed char __v16qs __attribute__ ((__vector_size__ (16)));
typedef unsigned char __v16qu __attribute__ ((__vector_size__ (16)));



typedef long long __m128i __attribute__ ((__vector_size__ (16), __may_alias__));
typedef double __m128d __attribute__ ((__vector_size__ (16), __may_alias__));


typedef long long __m128i_u __attribute__ ((__vector_size__ (16), __may_alias__, __aligned__ (1)));
typedef double __m128d_u __attribute__ ((__vector_size__ (16), __may_alias__, __aligned__ (1)));






extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_set_sd (double __F)
{
  return __extension__ (__m128d){ __F, 0.0 };
}


extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_set1_pd (double __F)
{
  return __extension__ (__m128d){ __F, __F };
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_set_pd1 (double __F)
{
  return _mm_set1_pd (__F);
}


extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_set_pd (double __W, double __X)
{
  return __extension__ (__m128d){ __X, __W };
}


extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_setr_pd (double __W, double __X)
{
  return __extension__ (__m128d){ __W, __X };
}


extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_undefined_pd (void)
{
  __m128d __Y = __Y;
  return __Y;
}


extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_setzero_pd (void)
{
  return __extension__ (__m128d){ 0.0, 0.0 };
}


extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_move_sd (__m128d __A, __m128d __B)
{
  return (__m128d) __builtin_ia32_movsd ((__v2df)__A, (__v2df)__B);
}


extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_load_pd (double const *__P)
{
  return *(__m128d *)__P;
}


extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_loadu_pd (double const *__P)
{
  return *(__m128d_u *)__P;
}


extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_load1_pd (double const *__P)
{
  return _mm_set1_pd (*__P);
}


extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_load_sd (double const *__P)
{
  return _mm_set_sd (*__P);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_load_pd1 (double const *__P)
{
  return _mm_load1_pd (__P);
}


extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_loadr_pd (double const *__P)
{
  __m128d __tmp = _mm_load_pd (__P);
  return __builtin_ia32_shufpd (__tmp, __tmp, (((0) << 1) | (1)));
}


extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_store_pd (double *__P, __m128d __A)
{
  *(__m128d *)__P = __A;
}


extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_storeu_pd (double *__P, __m128d __A)
{
  *(__m128d_u *)__P = __A;
}


extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_store_sd (double *__P, __m128d __A)
{
  *__P = ((__v2df)__A)[0];
}

extern __inline double __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtsd_f64 (__m128d __A)
{
  return ((__v2df)__A)[0];
}

extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_storel_pd (double *__P, __m128d __A)
{
  _mm_store_sd (__P, __A);
}


extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_storeh_pd (double *__P, __m128d __A)
{
  *__P = ((__v2df)__A)[1];
}



extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_store1_pd (double *__P, __m128d __A)
{
  _mm_store_pd (__P, __builtin_ia32_shufpd (__A, __A, (((0) << 1) | (0))));
}

extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_store_pd1 (double *__P, __m128d __A)
{
  _mm_store1_pd (__P, __A);
}


extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_storer_pd (double *__P, __m128d __A)
{
  _mm_store_pd (__P, __builtin_ia32_shufpd (__A, __A, (((0) << 1) | (1))));
}

extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtsi128_si32 (__m128i __A)
{
  return __builtin_ia32_vec_ext_v4si ((__v4si)__A, 0);
}



extern __inline long long __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtsi128_si64 (__m128i __A)
{
  return ((__v2di)__A)[0];
}


extern __inline long long __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtsi128_si64x (__m128i __A)
{
  return ((__v2di)__A)[0];
}


extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_add_pd (__m128d __A, __m128d __B)
{
  return (__m128d) ((__v2df)__A + (__v2df)__B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_add_sd (__m128d __A, __m128d __B)
{
  return (__m128d)__builtin_ia32_addsd ((__v2df)__A, (__v2df)__B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_sub_pd (__m128d __A, __m128d __B)
{
  return (__m128d) ((__v2df)__A - (__v2df)__B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_sub_sd (__m128d __A, __m128d __B)
{
  return (__m128d)__builtin_ia32_subsd ((__v2df)__A, (__v2df)__B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_mul_pd (__m128d __A, __m128d __B)
{
  return (__m128d) ((__v2df)__A * (__v2df)__B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_mul_sd (__m128d __A, __m128d __B)
{
  return (__m128d)__builtin_ia32_mulsd ((__v2df)__A, (__v2df)__B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_div_pd (__m128d __A, __m128d __B)
{
  return (__m128d) ((__v2df)__A / (__v2df)__B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_div_sd (__m128d __A, __m128d __B)
{
  return (__m128d)__builtin_ia32_divsd ((__v2df)__A, (__v2df)__B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_sqrt_pd (__m128d __A)
{
  return (__m128d)__builtin_ia32_sqrtpd ((__v2df)__A);
}


extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_sqrt_sd (__m128d __A, __m128d __B)
{
  __v2df __tmp = __builtin_ia32_movsd ((__v2df)__A, (__v2df)__B);
  return (__m128d)__builtin_ia32_sqrtsd ((__v2df)__tmp);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_min_pd (__m128d __A, __m128d __B)
{
  return (__m128d)__builtin_ia32_minpd ((__v2df)__A, (__v2df)__B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_min_sd (__m128d __A, __m128d __B)
{
  return (__m128d)__builtin_ia32_minsd ((__v2df)__A, (__v2df)__B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_max_pd (__m128d __A, __m128d __B)
{
  return (__m128d)__builtin_ia32_maxpd ((__v2df)__A, (__v2df)__B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_max_sd (__m128d __A, __m128d __B)
{
  return (__m128d)__builtin_ia32_maxsd ((__v2df)__A, (__v2df)__B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_and_pd (__m128d __A, __m128d __B)
{
  return (__m128d)__builtin_ia32_andpd ((__v2df)__A, (__v2df)__B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_andnot_pd (__m128d __A, __m128d __B)
{
  return (__m128d)__builtin_ia32_andnpd ((__v2df)__A, (__v2df)__B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_or_pd (__m128d __A, __m128d __B)
{
  return (__m128d)__builtin_ia32_orpd ((__v2df)__A, (__v2df)__B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_xor_pd (__m128d __A, __m128d __B)
{
  return (__m128d)__builtin_ia32_xorpd ((__v2df)__A, (__v2df)__B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpeq_pd (__m128d __A, __m128d __B)
{
  return (__m128d)__builtin_ia32_cmpeqpd ((__v2df)__A, (__v2df)__B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmplt_pd (__m128d __A, __m128d __B)
{
  return (__m128d)__builtin_ia32_cmpltpd ((__v2df)__A, (__v2df)__B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmple_pd (__m128d __A, __m128d __B)
{
  return (__m128d)__builtin_ia32_cmplepd ((__v2df)__A, (__v2df)__B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpgt_pd (__m128d __A, __m128d __B)
{
  return (__m128d)__builtin_ia32_cmpgtpd ((__v2df)__A, (__v2df)__B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpge_pd (__m128d __A, __m128d __B)
{
  return (__m128d)__builtin_ia32_cmpgepd ((__v2df)__A, (__v2df)__B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpneq_pd (__m128d __A, __m128d __B)
{
  return (__m128d)__builtin_ia32_cmpneqpd ((__v2df)__A, (__v2df)__B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpnlt_pd (__m128d __A, __m128d __B)
{
  return (__m128d)__builtin_ia32_cmpnltpd ((__v2df)__A, (__v2df)__B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpnle_pd (__m128d __A, __m128d __B)
{
  return (__m128d)__builtin_ia32_cmpnlepd ((__v2df)__A, (__v2df)__B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpngt_pd (__m128d __A, __m128d __B)
{
  return (__m128d)__builtin_ia32_cmpngtpd ((__v2df)__A, (__v2df)__B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpnge_pd (__m128d __A, __m128d __B)
{
  return (__m128d)__builtin_ia32_cmpngepd ((__v2df)__A, (__v2df)__B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpord_pd (__m128d __A, __m128d __B)
{
  return (__m128d)__builtin_ia32_cmpordpd ((__v2df)__A, (__v2df)__B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpunord_pd (__m128d __A, __m128d __B)
{
  return (__m128d)__builtin_ia32_cmpunordpd ((__v2df)__A, (__v2df)__B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpeq_sd (__m128d __A, __m128d __B)
{
  return (__m128d)__builtin_ia32_cmpeqsd ((__v2df)__A, (__v2df)__B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmplt_sd (__m128d __A, __m128d __B)
{
  return (__m128d)__builtin_ia32_cmpltsd ((__v2df)__A, (__v2df)__B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmple_sd (__m128d __A, __m128d __B)
{
  return (__m128d)__builtin_ia32_cmplesd ((__v2df)__A, (__v2df)__B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpgt_sd (__m128d __A, __m128d __B)
{
  return (__m128d) __builtin_ia32_movsd ((__v2df) __A,
      (__v2df)
      __builtin_ia32_cmpltsd ((__v2df) __B,
         (__v2df)
         __A));
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpge_sd (__m128d __A, __m128d __B)
{
  return (__m128d) __builtin_ia32_movsd ((__v2df) __A,
      (__v2df)
      __builtin_ia32_cmplesd ((__v2df) __B,
         (__v2df)
         __A));
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpneq_sd (__m128d __A, __m128d __B)
{
  return (__m128d)__builtin_ia32_cmpneqsd ((__v2df)__A, (__v2df)__B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpnlt_sd (__m128d __A, __m128d __B)
{
  return (__m128d)__builtin_ia32_cmpnltsd ((__v2df)__A, (__v2df)__B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpnle_sd (__m128d __A, __m128d __B)
{
  return (__m128d)__builtin_ia32_cmpnlesd ((__v2df)__A, (__v2df)__B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpngt_sd (__m128d __A, __m128d __B)
{
  return (__m128d) __builtin_ia32_movsd ((__v2df) __A,
      (__v2df)
      __builtin_ia32_cmpnltsd ((__v2df) __B,
          (__v2df)
          __A));
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpnge_sd (__m128d __A, __m128d __B)
{
  return (__m128d) __builtin_ia32_movsd ((__v2df) __A,
      (__v2df)
      __builtin_ia32_cmpnlesd ((__v2df) __B,
          (__v2df)
          __A));
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpord_sd (__m128d __A, __m128d __B)
{
  return (__m128d)__builtin_ia32_cmpordsd ((__v2df)__A, (__v2df)__B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpunord_sd (__m128d __A, __m128d __B)
{
  return (__m128d)__builtin_ia32_cmpunordsd ((__v2df)__A, (__v2df)__B);
}

extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_comieq_sd (__m128d __A, __m128d __B)
{
  return __builtin_ia32_comisdeq ((__v2df)__A, (__v2df)__B);
}

extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_comilt_sd (__m128d __A, __m128d __B)
{
  return __builtin_ia32_comisdlt ((__v2df)__A, (__v2df)__B);
}

extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_comile_sd (__m128d __A, __m128d __B)
{
  return __builtin_ia32_comisdle ((__v2df)__A, (__v2df)__B);
}

extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_comigt_sd (__m128d __A, __m128d __B)
{
  return __builtin_ia32_comisdgt ((__v2df)__A, (__v2df)__B);
}

extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_comige_sd (__m128d __A, __m128d __B)
{
  return __builtin_ia32_comisdge ((__v2df)__A, (__v2df)__B);
}

extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_comineq_sd (__m128d __A, __m128d __B)
{
  return __builtin_ia32_comisdneq ((__v2df)__A, (__v2df)__B);
}

extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_ucomieq_sd (__m128d __A, __m128d __B)
{
  return __builtin_ia32_ucomisdeq ((__v2df)__A, (__v2df)__B);
}

extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_ucomilt_sd (__m128d __A, __m128d __B)
{
  return __builtin_ia32_ucomisdlt ((__v2df)__A, (__v2df)__B);
}

extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_ucomile_sd (__m128d __A, __m128d __B)
{
  return __builtin_ia32_ucomisdle ((__v2df)__A, (__v2df)__B);
}

extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_ucomigt_sd (__m128d __A, __m128d __B)
{
  return __builtin_ia32_ucomisdgt ((__v2df)__A, (__v2df)__B);
}

extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_ucomige_sd (__m128d __A, __m128d __B)
{
  return __builtin_ia32_ucomisdge ((__v2df)__A, (__v2df)__B);
}

extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_ucomineq_sd (__m128d __A, __m128d __B)
{
  return __builtin_ia32_ucomisdneq ((__v2df)__A, (__v2df)__B);
}



extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_set_epi64x (long long __q1, long long __q0)
{
  return __extension__ (__m128i)(__v2di){ __q0, __q1 };
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_set_epi64 (__m64 __q1, __m64 __q0)
{
  return _mm_set_epi64x ((long long)__q1, (long long)__q0);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_set_epi32 (int __q3, int __q2, int __q1, int __q0)
{
  return __extension__ (__m128i)(__v4si){ __q0, __q1, __q2, __q3 };
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_set_epi16 (short __q7, short __q6, short __q5, short __q4,
        short __q3, short __q2, short __q1, short __q0)
{
  return __extension__ (__m128i)(__v8hi){
    __q0, __q1, __q2, __q3, __q4, __q5, __q6, __q7 };
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_set_epi8 (char __q15, char __q14, char __q13, char __q12,
       char __q11, char __q10, char __q09, char __q08,
       char __q07, char __q06, char __q05, char __q04,
       char __q03, char __q02, char __q01, char __q00)
{
  return __extension__ (__m128i)(__v16qi){
    __q00, __q01, __q02, __q03, __q04, __q05, __q06, __q07,
    __q08, __q09, __q10, __q11, __q12, __q13, __q14, __q15
  };
}



extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_set1_epi64x (long long __A)
{
  return _mm_set_epi64x (__A, __A);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_set1_epi64 (__m64 __A)
{
  return _mm_set_epi64 (__A, __A);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_set1_epi32 (int __A)
{
  return _mm_set_epi32 (__A, __A, __A, __A);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_set1_epi16 (short __A)
{
  return _mm_set_epi16 (__A, __A, __A, __A, __A, __A, __A, __A);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_set1_epi8 (char __A)
{
  return _mm_set_epi8 (__A, __A, __A, __A, __A, __A, __A, __A,
         __A, __A, __A, __A, __A, __A, __A, __A);
}




extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_setr_epi64 (__m64 __q0, __m64 __q1)
{
  return _mm_set_epi64 (__q1, __q0);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_setr_epi32 (int __q0, int __q1, int __q2, int __q3)
{
  return _mm_set_epi32 (__q3, __q2, __q1, __q0);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_setr_epi16 (short __q0, short __q1, short __q2, short __q3,
         short __q4, short __q5, short __q6, short __q7)
{
  return _mm_set_epi16 (__q7, __q6, __q5, __q4, __q3, __q2, __q1, __q0);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_setr_epi8 (char __q00, char __q01, char __q02, char __q03,
        char __q04, char __q05, char __q06, char __q07,
        char __q08, char __q09, char __q10, char __q11,
        char __q12, char __q13, char __q14, char __q15)
{
  return _mm_set_epi8 (__q15, __q14, __q13, __q12, __q11, __q10, __q09, __q08,
         __q07, __q06, __q05, __q04, __q03, __q02, __q01, __q00);
}



extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_load_si128 (__m128i const *__P)
{
  return *__P;
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_loadu_si128 (__m128i_u const *__P)
{
  return *__P;
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_loadl_epi64 (__m128i_u const *__P)
{
  return _mm_set_epi64 ((__m64)0LL, *(__m64_u *)__P);
}

extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_store_si128 (__m128i *__P, __m128i __B)
{
  *__P = __B;
}

extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_storeu_si128 (__m128i_u *__P, __m128i __B)
{
  *__P = __B;
}

extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_storel_epi64 (__m128i_u *__P, __m128i __B)
{
  *(__m64_u *)__P = (__m64) ((__v2di)__B)[0];
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_movepi64_pi64 (__m128i __B)
{
  return (__m64) ((__v2di)__B)[0];
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_movpi64_epi64 (__m64 __A)
{
  return _mm_set_epi64 ((__m64)0LL, __A);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_move_epi64 (__m128i __A)
{
  return (__m128i)__builtin_ia32_movq128 ((__v2di) __A);
}


extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_undefined_si128 (void)
{
  __m128i __Y = __Y;
  return __Y;
}


extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_setzero_si128 (void)
{
  return __extension__ (__m128i)(__v4si){ 0, 0, 0, 0 };
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtepi32_pd (__m128i __A)
{
  return (__m128d)__builtin_ia32_cvtdq2pd ((__v4si) __A);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtepi32_ps (__m128i __A)
{
  return (__m128)__builtin_ia32_cvtdq2ps ((__v4si) __A);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtpd_epi32 (__m128d __A)
{
  return (__m128i)__builtin_ia32_cvtpd2dq ((__v2df) __A);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtpd_pi32 (__m128d __A)
{
  return (__m64)__builtin_ia32_cvtpd2pi ((__v2df) __A);
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtpd_ps (__m128d __A)
{
  return (__m128)__builtin_ia32_cvtpd2ps ((__v2df) __A);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvttpd_epi32 (__m128d __A)
{
  return (__m128i)__builtin_ia32_cvttpd2dq ((__v2df) __A);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvttpd_pi32 (__m128d __A)
{
  return (__m64)__builtin_ia32_cvttpd2pi ((__v2df) __A);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtpi32_pd (__m64 __A)
{
  return (__m128d)__builtin_ia32_cvtpi2pd ((__v2si) __A);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtps_epi32 (__m128 __A)
{
  return (__m128i)__builtin_ia32_cvtps2dq ((__v4sf) __A);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvttps_epi32 (__m128 __A)
{
  return (__m128i)__builtin_ia32_cvttps2dq ((__v4sf) __A);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtps_pd (__m128 __A)
{
  return (__m128d)__builtin_ia32_cvtps2pd ((__v4sf) __A);
}

extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtsd_si32 (__m128d __A)
{
  return __builtin_ia32_cvtsd2si ((__v2df) __A);
}



extern __inline long long __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtsd_si64 (__m128d __A)
{
  return __builtin_ia32_cvtsd2si64 ((__v2df) __A);
}


extern __inline long long __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtsd_si64x (__m128d __A)
{
  return __builtin_ia32_cvtsd2si64 ((__v2df) __A);
}


extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvttsd_si32 (__m128d __A)
{
  return __builtin_ia32_cvttsd2si ((__v2df) __A);
}



extern __inline long long __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvttsd_si64 (__m128d __A)
{
  return __builtin_ia32_cvttsd2si64 ((__v2df) __A);
}


extern __inline long long __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvttsd_si64x (__m128d __A)
{
  return __builtin_ia32_cvttsd2si64 ((__v2df) __A);
}


extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtsd_ss (__m128 __A, __m128d __B)
{
  return (__m128)__builtin_ia32_cvtsd2ss ((__v4sf) __A, (__v2df) __B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtsi32_sd (__m128d __A, int __B)
{
  return (__m128d)__builtin_ia32_cvtsi2sd ((__v2df) __A, __B);
}



extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtsi64_sd (__m128d __A, long long __B)
{
  return (__m128d)__builtin_ia32_cvtsi642sd ((__v2df) __A, __B);
}


extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtsi64x_sd (__m128d __A, long long __B)
{
  return (__m128d)__builtin_ia32_cvtsi642sd ((__v2df) __A, __B);
}


extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtss_sd (__m128d __A, __m128 __B)
{
  return (__m128d)__builtin_ia32_cvtss2sd ((__v2df) __A, (__v4sf)__B);
}
# 919 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/emmintrin.h" 3 4
extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_unpackhi_pd (__m128d __A, __m128d __B)
{
  return (__m128d)__builtin_ia32_unpckhpd ((__v2df)__A, (__v2df)__B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_unpacklo_pd (__m128d __A, __m128d __B)
{
  return (__m128d)__builtin_ia32_unpcklpd ((__v2df)__A, (__v2df)__B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_loadh_pd (__m128d __A, double const *__B)
{
  return (__m128d)__builtin_ia32_loadhpd ((__v2df)__A, __B);
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_loadl_pd (__m128d __A, double const *__B)
{
  return (__m128d)__builtin_ia32_loadlpd ((__v2df)__A, __B);
}

extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_movemask_pd (__m128d __A)
{
  return __builtin_ia32_movmskpd ((__v2df)__A);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_packs_epi16 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_packsswb128 ((__v8hi)__A, (__v8hi)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_packs_epi32 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_packssdw128 ((__v4si)__A, (__v4si)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_packus_epi16 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_packuswb128 ((__v8hi)__A, (__v8hi)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_unpackhi_epi8 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_punpckhbw128 ((__v16qi)__A, (__v16qi)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_unpackhi_epi16 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_punpckhwd128 ((__v8hi)__A, (__v8hi)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_unpackhi_epi32 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_punpckhdq128 ((__v4si)__A, (__v4si)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_unpackhi_epi64 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_punpckhqdq128 ((__v2di)__A, (__v2di)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_unpacklo_epi8 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_punpcklbw128 ((__v16qi)__A, (__v16qi)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_unpacklo_epi16 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_punpcklwd128 ((__v8hi)__A, (__v8hi)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_unpacklo_epi32 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_punpckldq128 ((__v4si)__A, (__v4si)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_unpacklo_epi64 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_punpcklqdq128 ((__v2di)__A, (__v2di)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_add_epi8 (__m128i __A, __m128i __B)
{
  return (__m128i) ((__v16qu)__A + (__v16qu)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_add_epi16 (__m128i __A, __m128i __B)
{
  return (__m128i) ((__v8hu)__A + (__v8hu)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_add_epi32 (__m128i __A, __m128i __B)
{
  return (__m128i) ((__v4su)__A + (__v4su)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_add_epi64 (__m128i __A, __m128i __B)
{
  return (__m128i) ((__v2du)__A + (__v2du)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_adds_epi8 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_paddsb128 ((__v16qi)__A, (__v16qi)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_adds_epi16 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_paddsw128 ((__v8hi)__A, (__v8hi)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_adds_epu8 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_paddusb128 ((__v16qi)__A, (__v16qi)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_adds_epu16 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_paddusw128 ((__v8hi)__A, (__v8hi)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_sub_epi8 (__m128i __A, __m128i __B)
{
  return (__m128i) ((__v16qu)__A - (__v16qu)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_sub_epi16 (__m128i __A, __m128i __B)
{
  return (__m128i) ((__v8hu)__A - (__v8hu)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_sub_epi32 (__m128i __A, __m128i __B)
{
  return (__m128i) ((__v4su)__A - (__v4su)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_sub_epi64 (__m128i __A, __m128i __B)
{
  return (__m128i) ((__v2du)__A - (__v2du)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_subs_epi8 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_psubsb128 ((__v16qi)__A, (__v16qi)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_subs_epi16 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_psubsw128 ((__v8hi)__A, (__v8hi)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_subs_epu8 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_psubusb128 ((__v16qi)__A, (__v16qi)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_subs_epu16 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_psubusw128 ((__v8hi)__A, (__v8hi)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_madd_epi16 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_pmaddwd128 ((__v8hi)__A, (__v8hi)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_mulhi_epi16 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_pmulhw128 ((__v8hi)__A, (__v8hi)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_mullo_epi16 (__m128i __A, __m128i __B)
{
  return (__m128i) ((__v8hu)__A * (__v8hu)__B);
}

extern __inline __m64 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_mul_su32 (__m64 __A, __m64 __B)
{
  return (__m64)__builtin_ia32_pmuludq ((__v2si)__A, (__v2si)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_mul_epu32 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_pmuludq128 ((__v4si)__A, (__v4si)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_slli_epi16 (__m128i __A, int __B)
{
  return (__m128i)__builtin_ia32_psllwi128 ((__v8hi)__A, __B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_slli_epi32 (__m128i __A, int __B)
{
  return (__m128i)__builtin_ia32_pslldi128 ((__v4si)__A, __B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_slli_epi64 (__m128i __A, int __B)
{
  return (__m128i)__builtin_ia32_psllqi128 ((__v2di)__A, __B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_srai_epi16 (__m128i __A, int __B)
{
  return (__m128i)__builtin_ia32_psrawi128 ((__v8hi)__A, __B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_srai_epi32 (__m128i __A, int __B)
{
  return (__m128i)__builtin_ia32_psradi128 ((__v4si)__A, __B);
}
# 1206 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/emmintrin.h" 3 4
extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_srli_epi16 (__m128i __A, int __B)
{
  return (__m128i)__builtin_ia32_psrlwi128 ((__v8hi)__A, __B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_srli_epi32 (__m128i __A, int __B)
{
  return (__m128i)__builtin_ia32_psrldi128 ((__v4si)__A, __B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_srli_epi64 (__m128i __A, int __B)
{
  return (__m128i)__builtin_ia32_psrlqi128 ((__v2di)__A, __B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_sll_epi16 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_psllw128((__v8hi)__A, (__v8hi)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_sll_epi32 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_pslld128((__v4si)__A, (__v4si)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_sll_epi64 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_psllq128((__v2di)__A, (__v2di)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_sra_epi16 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_psraw128 ((__v8hi)__A, (__v8hi)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_sra_epi32 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_psrad128 ((__v4si)__A, (__v4si)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_srl_epi16 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_psrlw128 ((__v8hi)__A, (__v8hi)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_srl_epi32 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_psrld128 ((__v4si)__A, (__v4si)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_srl_epi64 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_psrlq128 ((__v2di)__A, (__v2di)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_and_si128 (__m128i __A, __m128i __B)
{
  return (__m128i) ((__v2du)__A & (__v2du)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_andnot_si128 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_pandn128 ((__v2di)__A, (__v2di)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_or_si128 (__m128i __A, __m128i __B)
{
  return (__m128i) ((__v2du)__A | (__v2du)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_xor_si128 (__m128i __A, __m128i __B)
{
  return (__m128i) ((__v2du)__A ^ (__v2du)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpeq_epi8 (__m128i __A, __m128i __B)
{
  return (__m128i) ((__v16qs)__A == (__v16qs)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpeq_epi16 (__m128i __A, __m128i __B)
{
  return (__m128i) ((__v8hi)__A == (__v8hi)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpeq_epi32 (__m128i __A, __m128i __B)
{
  return (__m128i) ((__v4si)__A == (__v4si)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmplt_epi8 (__m128i __A, __m128i __B)
{
  return (__m128i) ((__v16qs)__A < (__v16qs)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmplt_epi16 (__m128i __A, __m128i __B)
{
  return (__m128i) ((__v8hi)__A < (__v8hi)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmplt_epi32 (__m128i __A, __m128i __B)
{
  return (__m128i) ((__v4si)__A < (__v4si)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpgt_epi8 (__m128i __A, __m128i __B)
{
  return (__m128i) ((__v16qs)__A > (__v16qs)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpgt_epi16 (__m128i __A, __m128i __B)
{
  return (__m128i) ((__v8hi)__A > (__v8hi)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cmpgt_epi32 (__m128i __A, __m128i __B)
{
  return (__m128i) ((__v4si)__A > (__v4si)__B);
}
# 1370 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/emmintrin.h" 3 4
extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_max_epi16 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_pmaxsw128 ((__v8hi)__A, (__v8hi)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_max_epu8 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_pmaxub128 ((__v16qi)__A, (__v16qi)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_min_epi16 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_pminsw128 ((__v8hi)__A, (__v8hi)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_min_epu8 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_pminub128 ((__v16qi)__A, (__v16qi)__B);
}

extern __inline int __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_movemask_epi8 (__m128i __A)
{
  return __builtin_ia32_pmovmskb128 ((__v16qi)__A);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_mulhi_epu16 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_pmulhuw128 ((__v8hi)__A, (__v8hi)__B);
}
# 1433 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/emmintrin.h" 3 4
extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_maskmoveu_si128 (__m128i __A, __m128i __B, char *__C)
{
  __builtin_ia32_maskmovdqu ((__v16qi)__A, (__v16qi)__B, __C);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_avg_epu8 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_pavgb128 ((__v16qi)__A, (__v16qi)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_avg_epu16 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_pavgw128 ((__v8hi)__A, (__v8hi)__B);
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_sad_epu8 (__m128i __A, __m128i __B)
{
  return (__m128i)__builtin_ia32_psadbw128 ((__v16qi)__A, (__v16qi)__B);
}

extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_stream_si32 (int *__A, int __B)
{
  __builtin_ia32_movnti (__A, __B);
}


extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_stream_si64 (long long int *__A, long long int __B)
{
  __builtin_ia32_movnti64 (__A, __B);
}


extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_stream_si128 (__m128i *__A, __m128i __B)
{
  __builtin_ia32_movntdq ((__v2di *)__A, (__v2di)__B);
}

extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_stream_pd (double *__A, __m128d __B)
{
  __builtin_ia32_movntpd (__A, (__v2df)__B);
}

extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_clflush (void const *__A)
{
  __builtin_ia32_clflush (__A);
}

extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_lfence (void)
{
  __builtin_ia32_lfence ();
}

extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_mfence (void)
{
  __builtin_ia32_mfence ();
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtsi32_si128 (int __A)
{
  return _mm_set_epi32 (0, 0, 0, __A);
}



extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtsi64_si128 (long long __A)
{
  return _mm_set_epi64x (0, __A);
}


extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_cvtsi64x_si128 (long long __A)
{
  return _mm_set_epi64x (0, __A);
}




extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_castpd_ps(__m128d __A)
{
  return (__m128) __A;
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_castpd_si128(__m128d __A)
{
  return (__m128i) __A;
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_castps_pd(__m128 __A)
{
  return (__m128d) __A;
}

extern __inline __m128i __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_castps_si128(__m128 __A)
{
  return (__m128i) __A;
}

extern __inline __m128 __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_castsi128_ps(__m128i __A)
{
  return (__m128) __A;
}

extern __inline __m128d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_castsi128_pd(__m128i __A)
{
  return (__m128d) __A;
}
# 1253 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/xmmintrin.h" 2 3 4
# 1264 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/xmmintrin.h" 3 4
extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm_pause (void)
{
  __builtin_ia32_pause ();
}
# 1877 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h" 2
# 1916 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"

# 1916 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
static inline PetscErrorCode PetscMemcpy(void *a,const void *b,size_t n)
{







  ;

  if (a != b && n > 0) {
# 1952 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
    memcpy((char*)(a),(char*)(b),n);

  }
  return(0);
}
# 1980 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
static inline PetscErrorCode PetscMemzero(void *a,size_t n)
{
  if (n > 0) {
# 2001 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
      memset((char*)a,0,n);




  }
  return 0;
}
# 2210 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
extern PetscErrorCode MPIU_File_write_all(MPI_File,void*,PetscMPIInt,MPI_Datatype,MPI_Status*);
extern PetscErrorCode MPIU_File_read_all(MPI_File,void*,PetscMPIInt,MPI_Datatype,MPI_Status*);
# 2245 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
static inline PetscErrorCode PetscBLASIntCast(PetscInt a,PetscBLASInt *b)
{
  ;
  *b = (PetscBLASInt)(a);



  return(0);
}
# 2273 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
static inline PetscErrorCode PetscMPIIntCast(PetscInt a,PetscMPIInt *b)
{
  ;
  *b = (PetscMPIInt)(a);



  return(0);
}
# 2295 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
# 1 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/limits.h" 1 3 4
# 34 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/limits.h" 3 4
# 1 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/syslimits.h" 1 3 4






# 1 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/limits.h" 1 3 4
# 194 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/limits.h" 3 4
# 1 "/usr/include/limits.h" 1 3 4
# 144 "/usr/include/limits.h" 3 4
# 1 "/usr/include/bits/posix1_lim.h" 1 3 4
# 160 "/usr/include/bits/posix1_lim.h" 3 4
# 1 "/usr/include/bits/local_lim.h" 1 3 4
# 38 "/usr/include/bits/local_lim.h" 3 4
# 1 "/usr/include/linux/limits.h" 1 3 4
# 39 "/usr/include/bits/local_lim.h" 2 3 4
# 161 "/usr/include/bits/posix1_lim.h" 2 3 4
# 145 "/usr/include/limits.h" 2 3 4



# 1 "/usr/include/bits/posix2_lim.h" 1 3 4
# 149 "/usr/include/limits.h" 2 3 4
# 195 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/limits.h" 2 3 4
# 8 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/syslimits.h" 2 3 4
# 35 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/limits.h" 2 3 4
# 2296 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h" 2


# 1 "/usr/include/sys/param.h" 1 3 4
# 23 "/usr/include/sys/param.h" 3 4
# 1 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/stddef.h" 1 3 4
# 24 "/usr/include/sys/param.h" 2 3 4


# 1 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/limits.h" 1 3 4
# 27 "/usr/include/sys/param.h" 2 3 4

# 1 "/usr/include/signal.h" 1 3 4
# 30 "/usr/include/signal.h" 3 4


# 1 "/usr/include/bits/sigset.h" 1 3 4
# 103 "/usr/include/bits/sigset.h" 3 4

# 103 "/usr/include/bits/sigset.h" 3 4
extern int __sigismember (const __sigset_t *, int);
extern int __sigaddset (__sigset_t *, int);
extern int __sigdelset (__sigset_t *, int);
# 33 "/usr/include/signal.h" 2 3 4







typedef __sig_atomic_t sig_atomic_t;

# 57 "/usr/include/signal.h" 3 4
# 1 "/usr/include/bits/signum.h" 1 3 4
# 58 "/usr/include/signal.h" 2 3 4
# 75 "/usr/include/signal.h" 3 4
# 1 "/usr/include/time.h" 1 3 4
# 76 "/usr/include/signal.h" 2 3 4




# 1 "/usr/include/bits/siginfo.h" 1 3 4
# 24 "/usr/include/bits/siginfo.h" 3 4
# 1 "/usr/include/bits/wordsize.h" 1 3 4
# 25 "/usr/include/bits/siginfo.h" 2 3 4







typedef union sigval
  {
    int sival_int;
    void *sival_ptr;
  } sigval_t;
# 58 "/usr/include/bits/siginfo.h" 3 4
typedef __clock_t __sigchld_clock_t;



typedef struct
  {
    int si_signo;
    int si_errno;

    int si_code;

    union
      {
 int _pad[((128 / sizeof (int)) - 4)];


 struct
   {
     __pid_t si_pid;
     __uid_t si_uid;
   } _kill;


 struct
   {
     int si_tid;
     int si_overrun;
     sigval_t si_sigval;
   } _timer;


 struct
   {
     __pid_t si_pid;
     __uid_t si_uid;
     sigval_t si_sigval;
   } _rt;


 struct
   {
     __pid_t si_pid;
     __uid_t si_uid;
     int si_status;
     __sigchld_clock_t si_utime;
     __sigchld_clock_t si_stime;
   } _sigchld;


 struct
   {
     void *si_addr;
   } _sigfault;


 struct
   {
     long int si_band;
     int si_fd;
   } _sigpoll;


 struct
   {
     void *_call_addr;
     int _syscall;
     unsigned int _arch;
   } _sigsys;
      } _sifields;
  } siginfo_t ;
# 151 "/usr/include/bits/siginfo.h" 3 4
enum
{
  SI_ASYNCNL = -60,

  SI_TKILL = -6,

  SI_SIGIO,

  SI_ASYNCIO,

  SI_MESGQ,

  SI_TIMER,

  SI_QUEUE,

  SI_USER,

  SI_KERNEL = 0x80

};



enum
{
  ILL_ILLOPC = 1,

  ILL_ILLOPN,

  ILL_ILLADR,

  ILL_ILLTRP,

  ILL_PRVOPC,

  ILL_PRVREG,

  ILL_COPROC,

  ILL_BADSTK

};


enum
{
  FPE_INTDIV = 1,

  FPE_INTOVF,

  FPE_FLTDIV,

  FPE_FLTOVF,

  FPE_FLTUND,

  FPE_FLTRES,

  FPE_FLTINV,

  FPE_FLTSUB

};


enum
{
  SEGV_MAPERR = 1,

  SEGV_ACCERR

};


enum
{
  BUS_ADRALN = 1,

  BUS_ADRERR,

  BUS_OBJERR

};


enum
{
  TRAP_BRKPT = 1,

  TRAP_TRACE

};


enum
{
  CLD_EXITED = 1,

  CLD_KILLED,

  CLD_DUMPED,

  CLD_TRAPPED,

  CLD_STOPPED,

  CLD_CONTINUED

};


enum
{
  POLL_IN = 1,

  POLL_OUT,

  POLL_MSG,

  POLL_ERR,

  POLL_PRI,

  POLL_HUP

};
# 301 "/usr/include/bits/siginfo.h" 3 4
typedef struct sigevent
  {
    sigval_t sigev_value;
    int sigev_signo;
    int sigev_notify;

    union
      {
 int _pad[((64 / sizeof (int)) - 4)];



 __pid_t _tid;

 struct
   {
     void (*_function) (sigval_t);
     pthread_attr_t *_attribute;
   } _sigev_thread;
      } _sigev_un;
  } sigevent_t;






enum
{
  SIGEV_SIGNAL = 0,

  SIGEV_NONE,

  SIGEV_THREAD,


  SIGEV_THREAD_ID = 4

};
# 81 "/usr/include/signal.h" 2 3 4




typedef void (*__sighandler_t) (int);




extern __sighandler_t __sysv_signal (int __sig, __sighandler_t __handler)
     __attribute__ ((__nothrow__ , __leaf__));
# 100 "/usr/include/signal.h" 3 4


extern __sighandler_t signal (int __sig, __sighandler_t __handler)
     __attribute__ ((__nothrow__ , __leaf__));
# 114 "/usr/include/signal.h" 3 4

# 127 "/usr/include/signal.h" 3 4
extern int kill (__pid_t __pid, int __sig) __attribute__ ((__nothrow__ , __leaf__));






extern int killpg (__pid_t __pgrp, int __sig) __attribute__ ((__nothrow__ , __leaf__));




extern int raise (int __sig) __attribute__ ((__nothrow__ , __leaf__));




extern __sighandler_t ssignal (int __sig, __sighandler_t __handler)
     __attribute__ ((__nothrow__ , __leaf__));
extern int gsignal (int __sig) __attribute__ ((__nothrow__ , __leaf__));




extern void psignal (int __sig, const char *__s);




extern void psiginfo (const siginfo_t *__pinfo, const char *__s);
# 169 "/usr/include/signal.h" 3 4
extern int __sigpause (int __sig_or_mask, int __is_sig);
# 197 "/usr/include/signal.h" 3 4
extern int sigblock (int __mask) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__deprecated__));


extern int sigsetmask (int __mask) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__deprecated__));


extern int siggetmask (void) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__deprecated__));
# 217 "/usr/include/signal.h" 3 4
typedef __sighandler_t sig_t;





extern int sigemptyset (sigset_t *__set) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));


extern int sigfillset (sigset_t *__set) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));


extern int sigaddset (sigset_t *__set, int __signo) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));


extern int sigdelset (sigset_t *__set, int __signo) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));


extern int sigismember (const sigset_t *__set, int __signo)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));
# 253 "/usr/include/signal.h" 3 4
# 1 "/usr/include/bits/sigaction.h" 1 3 4
# 24 "/usr/include/bits/sigaction.h" 3 4
struct sigaction
  {


    union
      {

 __sighandler_t sa_handler;

 void (*sa_sigaction) (int, siginfo_t *, void *);
      }
    __sigaction_handler;







    __sigset_t sa_mask;


    int sa_flags;


    void (*sa_restorer) (void);
  };
# 254 "/usr/include/signal.h" 2 3 4


extern int sigprocmask (int __how, const sigset_t *__restrict __set,
   sigset_t *__restrict __oset) __attribute__ ((__nothrow__ , __leaf__));






extern int sigsuspend (const sigset_t *__set) __attribute__ ((__nonnull__ (1)));


extern int sigaction (int __sig, const struct sigaction *__restrict __act,
        struct sigaction *__restrict __oact) __attribute__ ((__nothrow__ , __leaf__));


extern int sigpending (sigset_t *__set) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));






extern int sigwait (const sigset_t *__restrict __set, int *__restrict __sig)
     __attribute__ ((__nonnull__ (1, 2)));






extern int sigwaitinfo (const sigset_t *__restrict __set,
   siginfo_t *__restrict __info) __attribute__ ((__nonnull__ (1)));






extern int sigtimedwait (const sigset_t *__restrict __set,
    siginfo_t *__restrict __info,
    const struct timespec *__restrict __timeout)
     __attribute__ ((__nonnull__ (1)));



extern int sigqueue (__pid_t __pid, int __sig, const union sigval __val)
     __attribute__ ((__nothrow__ , __leaf__));
# 311 "/usr/include/signal.h" 3 4
extern const char *const _sys_siglist[65];
extern const char *const sys_siglist[65];


struct sigvec
  {
    __sighandler_t sv_handler;
    int sv_mask;

    int sv_flags;

  };
# 335 "/usr/include/signal.h" 3 4
extern int sigvec (int __sig, const struct sigvec *__vec,
     struct sigvec *__ovec) __attribute__ ((__nothrow__ , __leaf__));



# 1 "/usr/include/bits/sigcontext.h" 1 3 4
# 29 "/usr/include/bits/sigcontext.h" 3 4
struct _fpx_sw_bytes
{
  __uint32_t magic1;
  __uint32_t extended_size;
  __uint64_t xstate_bv;
  __uint32_t xstate_size;
  __uint32_t padding[7];
};

struct _fpreg
{
  unsigned short significand[4];
  unsigned short exponent;
};

struct _fpxreg
{
  unsigned short significand[4];
  unsigned short exponent;
  unsigned short padding[3];
};

struct _xmmreg
{
  __uint32_t element[4];
};
# 121 "/usr/include/bits/sigcontext.h" 3 4
struct _fpstate
{

  __uint16_t cwd;
  __uint16_t swd;
  __uint16_t ftw;
  __uint16_t fop;
  __uint64_t rip;
  __uint64_t rdp;
  __uint32_t mxcsr;
  __uint32_t mxcr_mask;
  struct _fpxreg _st[8];
  struct _xmmreg _xmm[16];
  __uint32_t padding[24];
};

struct sigcontext
{
  __uint64_t r8;
  __uint64_t r9;
  __uint64_t r10;
  __uint64_t r11;
  __uint64_t r12;
  __uint64_t r13;
  __uint64_t r14;
  __uint64_t r15;
  __uint64_t rdi;
  __uint64_t rsi;
  __uint64_t rbp;
  __uint64_t rbx;
  __uint64_t rdx;
  __uint64_t rax;
  __uint64_t rcx;
  __uint64_t rsp;
  __uint64_t rip;
  __uint64_t eflags;
  unsigned short cs;
  unsigned short gs;
  unsigned short fs;
  unsigned short __pad0;
  __uint64_t err;
  __uint64_t trapno;
  __uint64_t oldmask;
  __uint64_t cr2;
  __extension__ union
    {
      struct _fpstate * fpstate;
      __uint64_t __fpstate_word;
    };
  __uint64_t __reserved1 [8];
};



struct _xsave_hdr
{
  __uint64_t xstate_bv;
  __uint64_t reserved1[2];
  __uint64_t reserved2[5];
};

struct _ymmh_state
{
  __uint32_t ymmh_space[64];
};

struct _xstate
{
  struct _fpstate fpstate;
  struct _xsave_hdr xstate_hdr;
  struct _ymmh_state ymmh;
};
# 341 "/usr/include/signal.h" 2 3 4


extern int sigreturn (struct sigcontext *__scp) __attribute__ ((__nothrow__ , __leaf__));






# 1 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/stddef.h" 1 3 4
# 351 "/usr/include/signal.h" 2 3 4




extern int siginterrupt (int __sig, int __interrupt) __attribute__ ((__nothrow__ , __leaf__));

# 1 "/usr/include/bits/sigstack.h" 1 3 4
# 25 "/usr/include/bits/sigstack.h" 3 4
struct sigstack
  {
    void *ss_sp;
    int ss_onstack;
  };



enum
{
  SS_ONSTACK = 1,

  SS_DISABLE

};
# 49 "/usr/include/bits/sigstack.h" 3 4
typedef struct sigaltstack
  {
    void *ss_sp;
    int ss_flags;
    size_t ss_size;
  } stack_t;
# 358 "/usr/include/signal.h" 2 3 4


# 1 "/usr/include/sys/ucontext.h" 1 3 4
# 22 "/usr/include/sys/ucontext.h" 3 4
# 1 "/usr/include/signal.h" 1 3 4
# 23 "/usr/include/sys/ucontext.h" 2 3 4
# 31 "/usr/include/sys/ucontext.h" 3 4
__extension__ typedef long long int greg_t;





typedef greg_t gregset_t[23];
# 92 "/usr/include/sys/ucontext.h" 3 4
struct _libc_fpxreg
{
  unsigned short int significand[4];
  unsigned short int exponent;
  unsigned short int padding[3];
};

struct _libc_xmmreg
{
  __uint32_t element[4];
};

struct _libc_fpstate
{

  __uint16_t cwd;
  __uint16_t swd;
  __uint16_t ftw;
  __uint16_t fop;
  __uint64_t rip;
  __uint64_t rdp;
  __uint32_t mxcsr;
  __uint32_t mxcr_mask;
  struct _libc_fpxreg _st[8];
  struct _libc_xmmreg _xmm[16];
  __uint32_t padding[24];
};


typedef struct _libc_fpstate *fpregset_t;


typedef struct
  {
    gregset_t gregs;

    fpregset_t fpregs;
    __extension__ unsigned long long __reserved1 [8];
} mcontext_t;


typedef struct ucontext
  {
    unsigned long int uc_flags;
    struct ucontext *uc_link;
    stack_t uc_stack;
    mcontext_t uc_mcontext;
    __sigset_t uc_sigmask;
    struct _libc_fpstate __fpregs_mem;
  } ucontext_t;
# 361 "/usr/include/signal.h" 2 3 4





extern int sigstack (struct sigstack *__ss, struct sigstack *__oss)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__deprecated__));



extern int sigaltstack (const struct sigaltstack *__restrict __ss,
   struct sigaltstack *__restrict __oss) __attribute__ ((__nothrow__ , __leaf__));
# 396 "/usr/include/signal.h" 3 4
# 1 "/usr/include/bits/sigthread.h" 1 3 4
# 30 "/usr/include/bits/sigthread.h" 3 4
extern int pthread_sigmask (int __how,
       const __sigset_t *__restrict __newmask,
       __sigset_t *__restrict __oldmask)__attribute__ ((__nothrow__ , __leaf__));


extern int pthread_kill (pthread_t __threadid, int __signo) __attribute__ ((__nothrow__ , __leaf__));
# 397 "/usr/include/signal.h" 2 3 4






extern int __libc_current_sigrtmin (void) __attribute__ ((__nothrow__ , __leaf__));

extern int __libc_current_sigrtmax (void) __attribute__ ((__nothrow__ , __leaf__));




# 29 "/usr/include/sys/param.h" 2 3 4


# 1 "/usr/include/bits/param.h" 1 3 4
# 28 "/usr/include/bits/param.h" 3 4
# 1 "/usr/include/linux/param.h" 1 3 4



# 1 "/usr/include/asm/param.h" 1 3 4
# 1 "/usr/include/asm-generic/param.h" 1 3 4
# 1 "/usr/include/asm/param.h" 2 3 4
# 5 "/usr/include/linux/param.h" 2 3 4
# 29 "/usr/include/bits/param.h" 2 3 4
# 32 "/usr/include/sys/param.h" 2 3 4
# 2299 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h" 2
# 2392 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"

# 2392 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
extern PetscErrorCode PetscGetArchType(char[],size_t);
extern PetscErrorCode PetscGetHostName(char[],size_t);
extern PetscErrorCode PetscGetUserName(char[],size_t);
extern PetscErrorCode PetscGetProgramName(char[],size_t);
extern PetscErrorCode PetscSetProgramName(const char[]);
extern PetscErrorCode PetscGetDate(char[],size_t);
extern PetscErrorCode PetscGetVersion(char[], size_t);

extern PetscErrorCode PetscSortInt(PetscInt,PetscInt[]);
extern PetscErrorCode PetscSortRemoveDupsInt(PetscInt*,PetscInt[]);
extern PetscErrorCode PetscFindInt(PetscInt, PetscInt, const PetscInt[], PetscInt*);
extern PetscErrorCode PetscSortIntWithPermutation(PetscInt,const PetscInt[],PetscInt[]);
extern PetscErrorCode PetscSortStrWithPermutation(PetscInt,const char*[],PetscInt[]);
extern PetscErrorCode PetscSortIntWithArray(PetscInt,PetscInt[],PetscInt[]);
extern PetscErrorCode PetscSortIntWithArrayPair(PetscInt,PetscInt*,PetscInt*,PetscInt*);
extern PetscErrorCode PetscSortMPIInt(PetscInt,PetscMPIInt[]);
extern PetscErrorCode PetscSortRemoveDupsMPIInt(PetscInt*,PetscMPIInt[]);
extern PetscErrorCode PetscSortMPIIntWithArray(PetscMPIInt,PetscMPIInt[],PetscMPIInt[]);
extern PetscErrorCode PetscSortIntWithScalarArray(PetscInt,PetscInt[],PetscScalar[]);
extern PetscErrorCode PetscSortIntWithDataArray(PetscInt,PetscInt[],void*,size_t,void*);
extern PetscErrorCode PetscSortReal(PetscInt,PetscReal[]);
extern PetscErrorCode PetscSortRealWithPermutation(PetscInt,const PetscReal[],PetscInt[]);
extern PetscErrorCode PetscSortRemoveDupsReal(PetscInt*,PetscReal[]);
extern PetscErrorCode PetscSortSplit(PetscInt,PetscInt,PetscScalar[],PetscInt[]);
extern PetscErrorCode PetscSortSplitReal(PetscInt,PetscInt,PetscReal[],PetscInt[]);
extern PetscErrorCode PetscProcessTree(PetscInt,const PetscBool [],const PetscInt[],PetscInt*,PetscInt**,PetscInt**,PetscInt**,PetscInt**);
extern PetscErrorCode PetscMergeIntArrayPair(PetscInt,const PetscInt*,const PetscInt*,PetscInt,const PetscInt*,const PetscInt*,PetscInt*,PetscInt**,PetscInt**);
extern PetscErrorCode PetscMergeIntArray(PetscInt,const PetscInt*,PetscInt,const PetscInt*,PetscInt*,PetscInt**);

extern PetscErrorCode PetscSetDisplay(void);
extern PetscErrorCode PetscGetDisplay(char[],size_t);
# 2434 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
typedef const char* PetscRandomType;





extern PetscClassId PETSC_RANDOM_CLASSID;

extern PetscErrorCode PetscRandomInitializePackage(void);
# 2453 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
typedef struct _p_PetscRandom* PetscRandom;


extern PetscFunctionList PetscRandomList;

extern PetscErrorCode PetscRandomRegister(const char[],PetscErrorCode (*)(PetscRandom));
extern PetscErrorCode PetscRandomSetType(PetscRandom, PetscRandomType);
extern PetscErrorCode PetscRandomSetFromOptions(PetscRandom);
extern PetscErrorCode PetscRandomGetType(PetscRandom, PetscRandomType*);
static inline PetscErrorCode PetscRandomViewFromOptions(PetscRandom A,PetscObject obj,const char name[]) {return PetscObjectViewFromOptions((PetscObject)A,obj,name);}
extern PetscErrorCode PetscRandomView(PetscRandom,PetscViewer);

extern PetscErrorCode PetscRandomCreate(MPI_Comm,PetscRandom*);
extern PetscErrorCode PetscRandomGetValue(PetscRandom,PetscScalar*);
extern PetscErrorCode PetscRandomGetValueReal(PetscRandom,PetscReal*);
extern PetscErrorCode PetscRandomGetInterval(PetscRandom,PetscScalar*,PetscScalar*);
extern PetscErrorCode PetscRandomSetInterval(PetscRandom,PetscScalar,PetscScalar);
extern PetscErrorCode PetscRandomSetSeed(PetscRandom,unsigned long);
extern PetscErrorCode PetscRandomGetSeed(PetscRandom,unsigned long *);
extern PetscErrorCode PetscRandomSeed(PetscRandom);
extern PetscErrorCode PetscRandomDestroy(PetscRandom*);

extern PetscErrorCode PetscGetFullPath(const char[],char[],size_t);
extern PetscErrorCode PetscGetRelativePath(const char[],char[],size_t);
extern PetscErrorCode PetscGetWorkingDirectory(char[],size_t);
extern PetscErrorCode PetscGetRealPath(const char[],char[]);
extern PetscErrorCode PetscGetHomeDirectory(char[],size_t);
extern PetscErrorCode PetscTestFile(const char[],char,PetscBool *);
extern PetscErrorCode PetscTestDirectory(const char[],char,PetscBool *);

extern PetscErrorCode PetscBinaryRead(int,void*,PetscInt,PetscDataType);
extern PetscErrorCode PetscBinarySynchronizedRead(MPI_Comm,int,void*,PetscInt,PetscDataType);
extern PetscErrorCode PetscBinarySynchronizedWrite(MPI_Comm,int,void*,PetscInt,PetscDataType,PetscBool );
extern PetscErrorCode PetscBinaryWrite(int,void*,PetscInt,PetscDataType,PetscBool );
extern PetscErrorCode PetscBinaryOpen(const char[],PetscFileMode,int *);
extern PetscErrorCode PetscBinaryClose(int);
extern PetscErrorCode PetscSharedTmp(MPI_Comm,PetscBool *);
extern PetscErrorCode PetscSharedWorkingDirectory(MPI_Comm,PetscBool *);
extern PetscErrorCode PetscGetTmp(MPI_Comm,char[],size_t);
extern PetscErrorCode PetscFileRetrieve(MPI_Comm,const char[],char[],size_t,PetscBool *);
extern PetscErrorCode PetscLs(MPI_Comm,const char[],char[],size_t,PetscBool *);
extern PetscErrorCode PetscOpenSocket(const char[],int,int*);
# 2516 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
typedef enum {PETSC_BINARY_SEEK_SET = 0,PETSC_BINARY_SEEK_CUR = 1,PETSC_BINARY_SEEK_END = 2} PetscBinarySeekType;
extern PetscErrorCode PetscBinarySeek(int,off_t,PetscBinarySeekType,off_t*);
extern PetscErrorCode PetscBinarySynchronizedSeek(MPI_Comm,int,off_t,PetscBinarySeekType,off_t*);
extern PetscErrorCode PetscByteSwap(void *,PetscDataType,PetscInt);

extern PetscErrorCode PetscSetDebugTerminal(const char[]);
extern PetscErrorCode PetscSetDebugger(const char[],PetscBool );
extern PetscErrorCode PetscSetDefaultDebugger(void);
extern PetscErrorCode PetscSetDebuggerFromString(const char*);
extern PetscErrorCode PetscAttachDebugger(void);
extern PetscErrorCode PetscStopForDebugger(void);

extern PetscErrorCode PetscGatherNumberOfMessages(MPI_Comm,const PetscMPIInt[],const PetscMPIInt[],PetscMPIInt*);
extern PetscErrorCode PetscGatherMessageLengths(MPI_Comm,PetscMPIInt,PetscMPIInt,const PetscMPIInt[],PetscMPIInt**,PetscMPIInt**);
extern PetscErrorCode PetscGatherMessageLengths2(MPI_Comm,PetscMPIInt,PetscMPIInt,const PetscMPIInt[],const PetscMPIInt[],PetscMPIInt**,PetscMPIInt**,PetscMPIInt**);
extern PetscErrorCode PetscPostIrecvInt(MPI_Comm,PetscMPIInt,PetscMPIInt,const PetscMPIInt[],const PetscMPIInt[],PetscInt***,MPI_Request**);
extern PetscErrorCode PetscPostIrecvScalar(MPI_Comm,PetscMPIInt,PetscMPIInt,const PetscMPIInt[],const PetscMPIInt[],PetscScalar***,MPI_Request**);
extern PetscErrorCode PetscCommBuildTwoSided(MPI_Comm,PetscMPIInt,MPI_Datatype,PetscInt,const PetscMPIInt*,const void*,PetscInt*,PetscMPIInt**,void*) ;
# 2548 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
typedef enum {
  PETSC_BUILDTWOSIDED_NOTSET = -1,
  PETSC_BUILDTWOSIDED_ALLREDUCE = 0,
  PETSC_BUILDTWOSIDED_IBARRIER = 1

} PetscBuildTwoSidedType;
extern const char *const PetscBuildTwoSidedTypes[];
extern PetscErrorCode PetscCommBuildTwoSidedSetType(MPI_Comm,PetscBuildTwoSidedType);
extern PetscErrorCode PetscCommBuildTwoSidedGetType(MPI_Comm,PetscBuildTwoSidedType*);

extern PetscErrorCode PetscSSEIsEnabled(MPI_Comm,PetscBool *,PetscBool *);
# 2569 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
 typedef enum {NOT_SET_VALUES, INSERT_VALUES, ADD_VALUES, MAX_VALUES, INSERT_ALL_VALUES, ADD_ALL_VALUES, INSERT_BC_VALUES, ADD_BC_VALUES} InsertMode;
# 2603 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
extern MPI_Comm PetscObjectComm(PetscObject);

typedef enum {PETSC_SUBCOMM_GENERAL=0,PETSC_SUBCOMM_CONTIGUOUS=1,PETSC_SUBCOMM_INTERLACED=2} PetscSubcommType;
extern const char *const PetscSubcommTypes[];
# 2642 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
typedef struct _n_PetscSubcomm* PetscSubcomm;

struct _n_PetscSubcomm {
  MPI_Comm parent;
  MPI_Comm dupparent;
  MPI_Comm child;
  PetscMPIInt n;
  PetscMPIInt color;
  PetscMPIInt *subsize;
  PetscSubcommType type;
};

static inline MPI_Comm PetscSubcommChild(PetscSubcomm scomm) {return scomm->child;}
static inline MPI_Comm PetscSubcommContiguousParent(PetscSubcomm scomm) {return scomm->dupparent;}
extern PetscErrorCode PetscSubcommCreate(MPI_Comm,PetscSubcomm*);
extern PetscErrorCode PetscSubcommDestroy(PetscSubcomm*);
extern PetscErrorCode PetscSubcommSetNumber(PetscSubcomm,PetscInt);
extern PetscErrorCode PetscSubcommSetType(PetscSubcomm,PetscSubcommType);
extern PetscErrorCode PetscSubcommSetTypeGeneral(PetscSubcomm,PetscMPIInt,PetscMPIInt);
extern PetscErrorCode PetscSubcommView(PetscSubcomm,PetscViewer);
extern PetscErrorCode PetscSubcommSetFromOptions(PetscSubcomm);
# 2671 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
typedef struct _n_PetscSegBuffer *PetscSegBuffer;
extern PetscErrorCode PetscSegBufferCreate(size_t,size_t,PetscSegBuffer*);
extern PetscErrorCode PetscSegBufferDestroy(PetscSegBuffer*);
extern PetscErrorCode PetscSegBufferGet(PetscSegBuffer,size_t,void*);
extern PetscErrorCode PetscSegBufferExtractAlloc(PetscSegBuffer,void*);
extern PetscErrorCode PetscSegBufferExtractTo(PetscSegBuffer,void*);
extern PetscErrorCode PetscSegBufferExtractInPlace(PetscSegBuffer,void*);
extern PetscErrorCode PetscSegBufferGetSize(PetscSegBuffer,size_t*);
extern PetscErrorCode PetscSegBufferUnuse(PetscSegBuffer,size_t);




static inline PetscErrorCode PetscSegBufferGetInts(PetscSegBuffer seg,PetscInt count,PetscInt *restrict *slot) {return PetscSegBufferGet(seg,(size_t)count,(void**)slot);}

extern PetscSegBuffer PetscCitationsList;
# 2703 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h"
static inline PetscErrorCode PetscCitationsRegister(const char cit[],PetscBool *set)
{
  size_t len;
  char *vstring;
  PetscErrorCode ierr;

  ;
  if (set && *set) return(0);
  ierr = PetscStrlen(cit,&len);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),2711,__func__,"/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = PetscSegBufferGet(PetscCitationsList,len,&vstring);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),2712,__func__,"/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = PetscMemcpy(vstring,cit,len);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),2713,__func__,"/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsys.h",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  if (set) *set = PETSC_TRUE;
  return(0);
}

extern PetscErrorCode PetscURLShorten(const char[],char[],size_t);
extern PetscErrorCode PetscGoogleDriveAuthorize(MPI_Comm,char[],char[],size_t);
extern PetscErrorCode PetscGoogleDriveRefresh(MPI_Comm,const char[],char[],size_t);
extern PetscErrorCode PetscGoogleDriveUpload(MPI_Comm,const char[],const char []);

extern PetscErrorCode PetscBoxAuthorize(MPI_Comm,char[],char[],size_t);
extern PetscErrorCode PetscBoxRefresh(MPI_Comm,const char[],char[],char[],size_t);

extern PetscErrorCode PetscTextBelt(MPI_Comm,const char[],const char[],PetscBool*);

extern PetscErrorCode PetscPullJSONValue(const char[],const char[],char[],size_t,PetscBool*);
extern PetscErrorCode PetscPushJSONValue(char[],const char[],const char[],size_t);
# 8 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscis.h" 2
# 1 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsftypes.h" 1
# 17 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsftypes.h"
typedef struct _p_PetscSF* PetscSF;
# 28 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsftypes.h"
typedef struct {
  PetscInt rank;
  PetscInt index;
} PetscSFNode;
# 9 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscis.h" 2
# 1 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscistypes.h" 1
# 13 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscistypes.h"
typedef struct _p_IS* IS;
# 32 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscistypes.h"
typedef struct _p_ISLocalToGlobalMapping* ISLocalToGlobalMapping;
# 49 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscistypes.h"
typedef struct _n_ISColoring* ISColoring;
# 59 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscistypes.h"
typedef struct _n_PetscLayout* PetscLayout;
# 75 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscistypes.h"
typedef struct _p_PetscSection *PetscSection;
# 10 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscis.h" 2


extern PetscClassId IS_CLASSID;

extern PetscErrorCode ISInitializePackage(void);
# 23 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscis.h"
typedef const char* ISType;





extern PetscFunctionList ISList;
extern PetscErrorCode ISSetType(IS, ISType);
extern PetscErrorCode ISGetType(IS, ISType *);
extern PetscErrorCode ISRegister(const char[],PetscErrorCode (*)(IS));
extern PetscErrorCode ISCreate(MPI_Comm,IS*);




extern PetscErrorCode ISCreateGeneral(MPI_Comm,PetscInt,const PetscInt[],PetscCopyMode,IS *);
extern PetscErrorCode ISGeneralSetIndices(IS,PetscInt,const PetscInt[],PetscCopyMode);
extern PetscErrorCode ISCreateBlock(MPI_Comm,PetscInt,PetscInt,const PetscInt[],PetscCopyMode,IS *);
extern PetscErrorCode ISBlockSetIndices(IS,PetscInt,PetscInt,const PetscInt[],PetscCopyMode);
extern PetscErrorCode ISCreateStride(MPI_Comm,PetscInt,PetscInt,PetscInt,IS *);
extern PetscErrorCode ISStrideSetStride(IS,PetscInt,PetscInt,PetscInt);

extern PetscErrorCode ISDestroy(IS*);
extern PetscErrorCode ISSetPermutation(IS);
extern PetscErrorCode ISPermutation(IS,PetscBool *);
extern PetscErrorCode ISSetIdentity(IS);
extern PetscErrorCode ISIdentity(IS,PetscBool *);
extern PetscErrorCode ISContiguousLocal(IS,PetscInt,PetscInt,PetscInt*,PetscBool*);

extern PetscErrorCode ISGetIndices(IS,const PetscInt *[]);
extern PetscErrorCode ISRestoreIndices(IS,const PetscInt *[]);
extern PetscErrorCode ISGetTotalIndices(IS,const PetscInt *[]);
extern PetscErrorCode ISRestoreTotalIndices(IS,const PetscInt *[]);
extern PetscErrorCode ISGetNonlocalIndices(IS,const PetscInt *[]);
extern PetscErrorCode ISRestoreNonlocalIndices(IS,const PetscInt *[]);
extern PetscErrorCode ISGetNonlocalIS(IS, IS *is);
extern PetscErrorCode ISRestoreNonlocalIS(IS, IS *is);
extern PetscErrorCode ISGetSize(IS,PetscInt *);
extern PetscErrorCode ISGetLocalSize(IS,PetscInt *);
extern PetscErrorCode ISInvertPermutation(IS,PetscInt,IS*);
extern PetscErrorCode ISView(IS,PetscViewer);
static inline PetscErrorCode ISViewFromOptions(IS A,PetscObject obj,const char name[]) {return PetscObjectViewFromOptions((PetscObject)A,obj,name);}
extern PetscErrorCode ISLoad(IS,PetscViewer);
extern PetscErrorCode ISEqual(IS,IS,PetscBool *);
extern PetscErrorCode ISSort(IS);
extern PetscErrorCode ISSortRemoveDups(IS);
extern PetscErrorCode ISSorted(IS,PetscBool *);
extern PetscErrorCode ISDifference(IS,IS,IS*);
extern PetscErrorCode ISSum(IS,IS,IS*);
extern PetscErrorCode ISExpand(IS,IS,IS*);
extern PetscErrorCode ISGetMinMax(IS,PetscInt*,PetscInt*);

extern PetscErrorCode ISBlockGetIndices(IS,const PetscInt *[]);
extern PetscErrorCode ISBlockRestoreIndices(IS,const PetscInt *[]);
extern PetscErrorCode ISBlockGetLocalSize(IS,PetscInt *);
extern PetscErrorCode ISBlockGetSize(IS,PetscInt *);
extern PetscErrorCode ISGetBlockSize(IS,PetscInt*);
extern PetscErrorCode ISSetBlockSize(IS,PetscInt);

extern PetscErrorCode ISStrideGetInfo(IS,PetscInt *,PetscInt*);

extern PetscErrorCode ISToGeneral(IS);

extern PetscErrorCode ISDuplicate(IS,IS*);
extern PetscErrorCode ISCopy(IS,IS);
extern PetscErrorCode ISAllGather(IS,IS*);
extern PetscErrorCode ISComplement(IS,PetscInt,PetscInt,IS*);
extern PetscErrorCode ISConcatenate(MPI_Comm,PetscInt,const IS[],IS*);
extern PetscErrorCode ISListToPair(MPI_Comm,PetscInt, IS[],IS*,IS*);
extern PetscErrorCode ISPairToList(IS,IS,PetscInt*, IS *[]);
extern PetscErrorCode ISEmbed(IS,IS,PetscBool,IS*);
extern PetscErrorCode ISSortPermutation(IS,PetscBool,IS*);
extern PetscErrorCode ISOnComm(IS,MPI_Comm,PetscCopyMode,IS*);


extern PetscClassId IS_LTOGM_CLASSID;
# 111 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscis.h"
typedef enum {IS_GTOLM_MASK,IS_GTOLM_DROP} ISGlobalToLocalMappingType;

extern PetscErrorCode ISLocalToGlobalMappingCreate(MPI_Comm,PetscInt,PetscInt,const PetscInt[],PetscCopyMode,ISLocalToGlobalMapping*);
extern PetscErrorCode ISLocalToGlobalMappingCreateIS(IS,ISLocalToGlobalMapping *);
extern PetscErrorCode ISLocalToGlobalMappingCreateSF(PetscSF,PetscInt,ISLocalToGlobalMapping*);
extern PetscErrorCode ISLocalToGlobalMappingView(ISLocalToGlobalMapping,PetscViewer);
extern PetscErrorCode ISLocalToGlobalMappingDestroy(ISLocalToGlobalMapping*);
extern PetscErrorCode ISLocalToGlobalMappingApply(ISLocalToGlobalMapping,PetscInt,const PetscInt[],PetscInt[]);
extern PetscErrorCode ISLocalToGlobalMappingApplyBlock(ISLocalToGlobalMapping,PetscInt,const PetscInt[],PetscInt[]);
extern PetscErrorCode ISLocalToGlobalMappingApplyIS(ISLocalToGlobalMapping,IS,IS*);
extern PetscErrorCode ISGlobalToLocalMappingApply(ISLocalToGlobalMapping,ISGlobalToLocalMappingType,PetscInt,const PetscInt[],PetscInt*,PetscInt[]);
extern PetscErrorCode ISGlobalToLocalMappingApplyBlock(ISLocalToGlobalMapping,ISGlobalToLocalMappingType,PetscInt,const PetscInt[],PetscInt*,PetscInt[]);
extern PetscErrorCode ISGlobalToLocalMappingApplyIS(ISLocalToGlobalMapping,ISGlobalToLocalMappingType,IS,IS*);
extern PetscErrorCode ISLocalToGlobalMappingGetSize(ISLocalToGlobalMapping,PetscInt*);
extern PetscErrorCode ISLocalToGlobalMappingGetInfo(ISLocalToGlobalMapping,PetscInt*,PetscInt*[],PetscInt*[],PetscInt**[]);
extern PetscErrorCode ISLocalToGlobalMappingRestoreInfo(ISLocalToGlobalMapping,PetscInt*,PetscInt*[],PetscInt*[],PetscInt**[]);
extern PetscErrorCode ISLocalToGlobalMappingGetBlockInfo(ISLocalToGlobalMapping,PetscInt*,PetscInt*[],PetscInt*[],PetscInt**[]);
extern PetscErrorCode ISLocalToGlobalMappingRestoreBlockInfo(ISLocalToGlobalMapping,PetscInt*,PetscInt*[],PetscInt*[],PetscInt**[]);
extern PetscErrorCode ISLocalToGlobalMappingGetIndices(ISLocalToGlobalMapping,const PetscInt**);
extern PetscErrorCode ISLocalToGlobalMappingRestoreIndices(ISLocalToGlobalMapping,const PetscInt**);
extern PetscErrorCode ISLocalToGlobalMappingGetBlockIndices(ISLocalToGlobalMapping,const PetscInt**);
extern PetscErrorCode ISLocalToGlobalMappingRestoreBlockIndices(ISLocalToGlobalMapping,const PetscInt**);
extern PetscErrorCode ISLocalToGlobalMappingConcatenate(MPI_Comm,PetscInt,const ISLocalToGlobalMapping[],ISLocalToGlobalMapping*);
extern PetscErrorCode ISG2LMapApply(ISLocalToGlobalMapping,PetscInt,const PetscInt[],PetscInt[]);
extern PetscErrorCode ISLocalToGlobalMappingGetBlockSize(ISLocalToGlobalMapping,PetscInt*);
# 155 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscis.h"
typedef enum {IS_COLORING_GLOBAL,IS_COLORING_GHOSTED} ISColoringType;
extern const char *const ISColoringTypes[];
typedef unsigned short ISColoringValue;
extern PetscErrorCode ISAllGatherColors(MPI_Comm,PetscInt,ISColoringValue*,PetscInt*,ISColoringValue*[]);

extern PetscErrorCode ISColoringCreate(MPI_Comm,PetscInt,PetscInt,const ISColoringValue[],PetscCopyMode,ISColoring*);
extern PetscErrorCode ISColoringDestroy(ISColoring*);
extern PetscErrorCode ISColoringView(ISColoring,PetscViewer);
extern PetscErrorCode ISColoringViewFromOptions(ISColoring,PetscObject,const char[]);
extern PetscErrorCode ISColoringGetIS(ISColoring,PetscInt*,IS*[]);
extern PetscErrorCode ISColoringRestoreIS(ISColoring,IS*[]);
extern PetscErrorCode ISColoringReference(ISColoring);
extern PetscErrorCode ISColoringSetType(ISColoring,ISColoringType);




extern PetscErrorCode ISPartitioningToNumbering(IS,IS*);
extern PetscErrorCode ISPartitioningCount(IS,PetscInt,PetscInt[]);

extern PetscErrorCode ISCompressIndicesGeneral(PetscInt,PetscInt,PetscInt,PetscInt,const IS[],IS[]);
extern PetscErrorCode ISCompressIndicesSorted(PetscInt,PetscInt,PetscInt,const IS[],IS[]);
extern PetscErrorCode ISExpandIndicesGeneral(PetscInt,PetscInt,PetscInt,PetscInt,const IS[],IS[]);


struct _n_PetscLayout{
  MPI_Comm comm;
  PetscInt n,N;
  PetscInt rstart,rend;
  PetscInt *range;
  PetscInt bs;



  PetscInt refcnt;
  ISLocalToGlobalMapping mapping;
  PetscInt *trstarts;
};
# 214 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscis.h"
static inline PetscErrorCode PetscLayoutFindOwner(PetscLayout map,PetscInt idx,PetscInt *owner)
{
  PetscErrorCode ierr;
  PetscMPIInt lo = 0,hi,t;

  ;
  *owner = -1;
  if (!((map->n >= 0) && (map->N >= 0) && (map->range))) return PetscError(((MPI_Comm)0x44000001),221,__func__,"/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscis.h",73,PETSC_ERROR_INITIAL,"PetscLayoutSetUp() must be called first");
  if (idx < 0 || idx > map->N) return PetscError(((MPI_Comm)0x44000001),222,__func__,"/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscis.h",63,PETSC_ERROR_INITIAL,"Index %D is out of range",idx);
  ierr = MPI_Comm_size(map->comm,&hi);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),223,__func__,"/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscis.h",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  while (hi - lo > 1) {
    t = lo + (hi - lo) / 2;
    if (idx < map->range[t]) hi = t;
    else lo = t;
  }
  *owner = lo;
  return(0);
}
# 254 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscis.h"
static inline PetscErrorCode PetscLayoutFindOwnerIndex(PetscLayout map,PetscInt idx,PetscInt *owner, PetscInt *lidx)
{
  PetscErrorCode ierr;
  PetscMPIInt lo = 0,hi,t;

  ;
  if (!((map->n >= 0) && (map->N >= 0) && (map->range))) return PetscError(((MPI_Comm)0x44000001),260,__func__,"/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscis.h",73,PETSC_ERROR_INITIAL,"PetscLayoutSetUp() must be called first");
  if (idx < 0 || idx > map->N) return PetscError(((MPI_Comm)0x44000001),261,__func__,"/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscis.h",63,PETSC_ERROR_INITIAL,"Index %D is out of range",idx);
  ierr = MPI_Comm_size(map->comm,&hi);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),262,__func__,"/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscis.h",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  while (hi - lo > 1) {
    t = lo + (hi - lo) / 2;
    if (idx < map->range[t]) hi = t;
    else lo = t;
  }
  if (owner) *owner = lo;
  if (lidx) *lidx = idx-map->range[lo];
  return(0);
}

extern PetscErrorCode PetscLayoutCreate(MPI_Comm,PetscLayout*);
extern PetscErrorCode PetscLayoutSetUp(PetscLayout);
extern PetscErrorCode PetscLayoutDestroy(PetscLayout*);
extern PetscErrorCode PetscLayoutDuplicate(PetscLayout,PetscLayout*);
extern PetscErrorCode PetscLayoutReference(PetscLayout,PetscLayout*);
extern PetscErrorCode PetscLayoutSetLocalSize(PetscLayout,PetscInt);
extern PetscErrorCode PetscLayoutGetLocalSize(PetscLayout,PetscInt *);
extern PetscErrorCode PetscLayoutSetSize(PetscLayout,PetscInt);
extern PetscErrorCode PetscLayoutGetSize(PetscLayout,PetscInt *);
extern PetscErrorCode PetscLayoutSetBlockSize(PetscLayout,PetscInt);
extern PetscErrorCode PetscLayoutGetBlockSize(PetscLayout,PetscInt*);
extern PetscErrorCode PetscLayoutGetRange(PetscLayout,PetscInt *,PetscInt *);
extern PetscErrorCode PetscLayoutGetRanges(PetscLayout,const PetscInt *[]);
extern PetscErrorCode PetscLayoutSetISLocalToGlobalMapping(PetscLayout,ISLocalToGlobalMapping);
extern PetscErrorCode PetscSFSetGraphLayout(PetscSF,PetscLayout,PetscInt,const PetscInt*,PetscCopyMode,const PetscInt*);

extern PetscClassId PETSC_SECTION_CLASSID;

extern PetscErrorCode PetscSectionCreate(MPI_Comm,PetscSection*);
extern PetscErrorCode PetscSectionClone(PetscSection, PetscSection*);
extern PetscErrorCode PetscSectionCopy(PetscSection, PetscSection);
extern PetscErrorCode PetscSectionGetNumFields(PetscSection, PetscInt *);
extern PetscErrorCode PetscSectionSetNumFields(PetscSection, PetscInt);
extern PetscErrorCode PetscSectionGetFieldName(PetscSection, PetscInt, const char *[]);
extern PetscErrorCode PetscSectionSetFieldName(PetscSection, PetscInt, const char []);
extern PetscErrorCode PetscSectionGetFieldComponents(PetscSection, PetscInt, PetscInt *);
extern PetscErrorCode PetscSectionSetFieldComponents(PetscSection, PetscInt, PetscInt);
extern PetscErrorCode PetscSectionGetChart(PetscSection, PetscInt *, PetscInt *);
extern PetscErrorCode PetscSectionSetChart(PetscSection, PetscInt, PetscInt);
extern PetscErrorCode PetscSectionGetPermutation(PetscSection, IS *);
extern PetscErrorCode PetscSectionSetPermutation(PetscSection, IS);
extern PetscErrorCode PetscSectionGetDof(PetscSection, PetscInt, PetscInt*);
extern PetscErrorCode PetscSectionSetDof(PetscSection, PetscInt, PetscInt);
extern PetscErrorCode PetscSectionAddDof(PetscSection, PetscInt, PetscInt);
extern PetscErrorCode PetscSectionGetFieldDof(PetscSection, PetscInt, PetscInt, PetscInt*);
extern PetscErrorCode PetscSectionSetFieldDof(PetscSection, PetscInt, PetscInt, PetscInt);
extern PetscErrorCode PetscSectionAddFieldDof(PetscSection, PetscInt, PetscInt, PetscInt);
extern PetscErrorCode PetscSectionHasConstraints(PetscSection, PetscBool *);
extern PetscErrorCode PetscSectionGetConstraintDof(PetscSection, PetscInt, PetscInt*);
extern PetscErrorCode PetscSectionSetConstraintDof(PetscSection, PetscInt, PetscInt);
extern PetscErrorCode PetscSectionAddConstraintDof(PetscSection, PetscInt, PetscInt);
extern PetscErrorCode PetscSectionGetFieldConstraintDof(PetscSection, PetscInt, PetscInt, PetscInt*);
extern PetscErrorCode PetscSectionSetFieldConstraintDof(PetscSection, PetscInt, PetscInt, PetscInt);
extern PetscErrorCode PetscSectionAddFieldConstraintDof(PetscSection, PetscInt, PetscInt, PetscInt);
extern PetscErrorCode PetscSectionGetConstraintIndices(PetscSection, PetscInt, const PetscInt**);
extern PetscErrorCode PetscSectionSetConstraintIndices(PetscSection, PetscInt, const PetscInt*);
extern PetscErrorCode PetscSectionGetFieldConstraintIndices(PetscSection, PetscInt, PetscInt, const PetscInt**);
extern PetscErrorCode PetscSectionSetFieldConstraintIndices(PetscSection, PetscInt, PetscInt, const PetscInt*);
extern PetscErrorCode PetscSectionSetUpBC(PetscSection);
extern PetscErrorCode PetscSectionSetUp(PetscSection);
extern PetscErrorCode PetscSectionGetMaxDof(PetscSection, PetscInt*);
extern PetscErrorCode PetscSectionGetStorageSize(PetscSection, PetscInt*);
extern PetscErrorCode PetscSectionGetConstrainedStorageSize(PetscSection, PetscInt*);
extern PetscErrorCode PetscSectionGetOffset(PetscSection, PetscInt, PetscInt*);
extern PetscErrorCode PetscSectionSetOffset(PetscSection, PetscInt, PetscInt);
extern PetscErrorCode PetscSectionGetFieldOffset(PetscSection, PetscInt, PetscInt, PetscInt*);
extern PetscErrorCode PetscSectionSetFieldOffset(PetscSection, PetscInt, PetscInt, PetscInt);
extern PetscErrorCode PetscSectionGetOffsetRange(PetscSection, PetscInt *, PetscInt *);
extern PetscErrorCode PetscSectionView(PetscSection, PetscViewer);
static inline PetscErrorCode PetscSectionViewFromOptions(PetscSection A,PetscObject obj,const char name[]) {return PetscObjectViewFromOptions((PetscObject)A,obj,name);}
extern PetscErrorCode PetscSectionReset(PetscSection);
extern PetscErrorCode PetscSectionDestroy(PetscSection*);
extern PetscErrorCode PetscSectionCreateGlobalSection(PetscSection, PetscSF, PetscBool, PetscSection *);
extern PetscErrorCode PetscSectionCreateGlobalSectionCensored(PetscSection, PetscSF, PetscBool, PetscInt, const PetscInt [], PetscSection *);
extern PetscErrorCode PetscSectionCreateSubsection(PetscSection, PetscInt, PetscInt [], PetscSection *);
extern PetscErrorCode PetscSectionCreateSubmeshSection(PetscSection, IS, PetscSection *);
extern PetscErrorCode PetscSectionGetPointLayout(MPI_Comm, PetscSection, PetscLayout *);
extern PetscErrorCode PetscSectionGetValueLayout(MPI_Comm, PetscSection, PetscLayout *);
extern PetscErrorCode PetscSectionPermute(PetscSection, IS, PetscSection *);
extern PetscErrorCode PetscSectionGetField(PetscSection, PetscInt, PetscSection *);

extern PetscErrorCode PetscSectionSetClosureIndex(PetscSection, PetscObject, PetscSection, IS);
extern PetscErrorCode PetscSectionGetClosureIndex(PetscSection, PetscObject, PetscSection *, IS *);


extern PetscErrorCode PetscSFConvertPartition(PetscSF, PetscSection, IS, ISLocalToGlobalMapping *, PetscSF *);
extern PetscErrorCode PetscSFCreateRemoteOffsets(PetscSF, PetscSection, PetscSection, PetscInt **);
extern PetscErrorCode PetscSFDistributeSection(PetscSF, PetscSection, PetscInt **, PetscSection);
extern PetscErrorCode PetscSFCreateSectionSF(PetscSF, PetscSection, PetscInt [], PetscSection, PetscSF *);
# 10 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscvec.h" 2
# 1 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscviewer.h" 1
# 11 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscviewer.h"
extern PetscClassId PETSC_VIEWER_CLASSID;
# 20 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscviewer.h"
typedef const char* PetscViewerType;
# 34 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscviewer.h"
extern PetscFunctionList PetscViewerList;
extern PetscErrorCode PetscViewerInitializePackage(void);

extern PetscErrorCode PetscViewerRegister(const char[],PetscErrorCode (*)(PetscViewer));

extern PetscErrorCode PetscViewerCreate(MPI_Comm,PetscViewer*);
extern PetscErrorCode PetscViewerSetFromOptions(PetscViewer);
extern PetscErrorCode PetscViewerASCIIOpenWithFILE(MPI_Comm,FILE*,PetscViewer*);

extern PetscErrorCode PetscViewerASCIIOpen(MPI_Comm,const char[],PetscViewer*);
extern PetscErrorCode PetscViewerASCIISetFILE(PetscViewer,FILE*);
extern PetscErrorCode PetscViewerBinaryOpen(MPI_Comm,const char[],PetscFileMode,PetscViewer*);
extern PetscErrorCode PetscViewerBinaryGetFlowControl(PetscViewer,PetscInt*);
extern PetscErrorCode PetscViewerBinarySetFlowControl(PetscViewer,PetscInt);
extern PetscErrorCode PetscViewerBinarySetUseMPIIO(PetscViewer,PetscBool);
extern PetscErrorCode PetscViewerBinaryGetUseMPIIO(PetscViewer,PetscBool *);

extern PetscErrorCode PetscViewerBinaryGetMPIIODescriptor(PetscViewer,MPI_File*);
extern PetscErrorCode PetscViewerBinaryGetMPIIOOffset(PetscViewer,MPI_Offset*);
extern PetscErrorCode PetscViewerBinaryAddMPIIOOffset(PetscViewer,MPI_Offset);


extern PetscErrorCode PetscViewerSocketOpen(MPI_Comm,const char[],int,PetscViewer*);
extern PetscErrorCode PetscViewerStringOpen(MPI_Comm,char[],size_t,PetscViewer*);
extern PetscErrorCode PetscViewerDrawOpen(MPI_Comm,const char[],const char[],int,int,int,int,PetscViewer*);
# 1 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscdrawtypes.h" 1
# 11 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscdrawtypes.h"
typedef const char* PetscDrawType;
# 28 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscdrawtypes.h"
typedef struct _p_PetscDraw* PetscDraw;
# 39 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscdrawtypes.h"
typedef struct _p_PetscDrawAxis* PetscDrawAxis;
# 50 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscdrawtypes.h"
typedef struct _p_PetscDrawLG* PetscDrawLG;
# 61 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscdrawtypes.h"
typedef struct _p_PetscDrawSP* PetscDrawSP;
# 72 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscdrawtypes.h"
typedef struct _p_PetscDrawHG* PetscDrawHG;
# 83 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscdrawtypes.h"
typedef struct _p_PetscDrawBar* PetscDrawBar;
# 60 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscviewer.h" 2
extern PetscErrorCode PetscViewerDrawSetDrawType(PetscViewer,PetscDrawType);
extern PetscErrorCode PetscViewerMathematicaOpen(MPI_Comm, int, const char[], const char[], PetscViewer *);
extern PetscErrorCode PetscViewerSiloOpen(MPI_Comm, const char[], PetscViewer *);
extern PetscErrorCode PetscViewerMatlabOpen(MPI_Comm,const char[],PetscFileMode,PetscViewer*);

extern PetscErrorCode PetscViewerGetType(PetscViewer,PetscViewerType*);
extern PetscErrorCode PetscViewerSetType(PetscViewer,PetscViewerType);
extern PetscErrorCode PetscViewerDestroy(PetscViewer*);
extern PetscErrorCode PetscViewerGetSingleton(PetscViewer,PetscViewer*);
extern PetscErrorCode PetscViewerRestoreSingleton(PetscViewer,PetscViewer*);
extern PetscErrorCode PetscViewerGetSubcomm(PetscViewer,MPI_Comm,PetscViewer*);
extern PetscErrorCode PetscViewerRestoreSubcomm(PetscViewer,MPI_Comm,PetscViewer*);

extern PetscErrorCode PetscViewerSetUp(PetscViewer);
extern PetscErrorCode PetscViewerView(PetscViewer,PetscViewer);
static inline PetscErrorCode PetscViewerViewFromOptions(PetscViewer A,PetscObject obj,const char name[]) {return PetscObjectViewFromOptions((PetscObject)A,obj,name);}

extern PetscErrorCode PetscViewerSetOptionsPrefix(PetscViewer,const char[]);
extern PetscErrorCode PetscViewerAppendOptionsPrefix(PetscViewer,const char[]);
extern PetscErrorCode PetscViewerGetOptionsPrefix(PetscViewer,const char*[]);
# 91 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscviewer.h"
typedef enum {
  PETSC_VIEWER_DEFAULT,
  PETSC_VIEWER_ASCII_MATLAB,
  PETSC_VIEWER_ASCII_MATHEMATICA,
  PETSC_VIEWER_ASCII_IMPL,
  PETSC_VIEWER_ASCII_INFO,
  PETSC_VIEWER_ASCII_INFO_DETAIL,
  PETSC_VIEWER_ASCII_COMMON,
  PETSC_VIEWER_ASCII_SYMMODU,
  PETSC_VIEWER_ASCII_INDEX,
  PETSC_VIEWER_ASCII_DENSE,
  PETSC_VIEWER_ASCII_MATRIXMARKET,
  PETSC_VIEWER_ASCII_VTK,
  PETSC_VIEWER_ASCII_VTK_CELL,
  PETSC_VIEWER_ASCII_VTK_COORDS,
  PETSC_VIEWER_ASCII_PCICE,
  PETSC_VIEWER_ASCII_PYTHON,
  PETSC_VIEWER_ASCII_FACTOR_INFO,
  PETSC_VIEWER_ASCII_LATEX,
  PETSC_VIEWER_DRAW_BASIC,
  PETSC_VIEWER_DRAW_LG,
  PETSC_VIEWER_DRAW_CONTOUR,
  PETSC_VIEWER_DRAW_PORTS,
  PETSC_VIEWER_VTK_VTS,
  PETSC_VIEWER_VTK_VTR,
  PETSC_VIEWER_VTK_VTU,
  PETSC_VIEWER_BINARY_MATLAB,
  PETSC_VIEWER_NATIVE,
  PETSC_VIEWER_HDF5_VIZ,
  PETSC_VIEWER_NOFORMAT
  } PetscViewerFormat;
extern const char *const PetscViewerFormats[];

extern PetscErrorCode PetscViewerSetFormat(PetscViewer,PetscViewerFormat);
extern PetscErrorCode PetscViewerPushFormat(PetscViewer,PetscViewerFormat);
extern PetscErrorCode PetscViewerPopFormat(PetscViewer);
extern PetscErrorCode PetscViewerGetFormat(PetscViewer,PetscViewerFormat*);
extern PetscErrorCode PetscViewerFlush(PetscViewer);

extern PetscErrorCode PetscOptionsGetViewer(MPI_Comm,const char[],const char[],PetscViewer*,PetscViewerFormat*,PetscBool*);

extern PetscErrorCode PetscOptionsViewer_Private(PetscOptions*,const char[],const char[],const char[],PetscViewer*,PetscViewerFormat *,PetscBool *);





extern PetscErrorCode PetscViewerASCIIGetPointer(PetscViewer,FILE**);
extern PetscErrorCode PetscViewerFileGetMode(PetscViewer,PetscFileMode*);
extern PetscErrorCode PetscViewerFileSetMode(PetscViewer,PetscFileMode);
extern PetscErrorCode PetscViewerRead(PetscViewer,void*,PetscInt,PetscInt*,PetscDataType);
extern PetscErrorCode PetscViewerASCIIPrintf(PetscViewer,const char[],...);
extern PetscErrorCode PetscViewerASCIISynchronizedPrintf(PetscViewer,const char[],...);
extern PetscErrorCode PetscViewerASCIISynchronizedAllow(PetscViewer,PetscBool);
extern PetscErrorCode PetscViewerASCIIPushTab(PetscViewer);
extern PetscErrorCode PetscViewerASCIIPopTab(PetscViewer);
extern PetscErrorCode PetscViewerASCIIUseTabs(PetscViewer,PetscBool );
extern PetscErrorCode PetscViewerASCIISetTab(PetscViewer,PetscInt);
extern PetscErrorCode PetscViewerASCIIGetTab(PetscViewer,PetscInt*);
extern PetscErrorCode PetscViewerASCIIAddTab(PetscViewer,PetscInt);
extern PetscErrorCode PetscViewerASCIISubtractTab(PetscViewer,PetscInt);
extern PetscErrorCode PetscViewerASCIIRead(PetscViewer,void *,PetscInt,PetscInt*,PetscDataType);
extern PetscErrorCode PetscViewerBinaryGetDescriptor(PetscViewer,int*);
extern PetscErrorCode PetscViewerBinaryGetInfoPointer(PetscViewer,FILE **);
extern PetscErrorCode PetscViewerBinaryRead(PetscViewer,void*,PetscInt,PetscInt*,PetscDataType);
extern PetscErrorCode PetscViewerBinaryWrite(PetscViewer,void*,PetscInt,PetscDataType,PetscBool );
extern PetscErrorCode PetscViewerStringSPrintf(PetscViewer,const char[],...);
extern PetscErrorCode PetscViewerStringSetString(PetscViewer,char[],PetscInt);
extern PetscErrorCode PetscViewerDrawClear(PetscViewer);
extern PetscErrorCode PetscViewerDrawSetHold(PetscViewer,PetscBool);
extern PetscErrorCode PetscViewerDrawGetHold(PetscViewer,PetscBool*);
extern PetscErrorCode PetscViewerDrawSetPause(PetscViewer,PetscReal);
extern PetscErrorCode PetscViewerDrawGetPause(PetscViewer,PetscReal*);
extern PetscErrorCode PetscViewerDrawSetInfo(PetscViewer,const char[],const char[],int,int,int,int);
extern PetscErrorCode PetscViewerDrawResize(PetscViewer,int,int);
extern PetscErrorCode PetscViewerDrawSetBounds(PetscViewer,PetscInt,const PetscReal*);
extern PetscErrorCode PetscViewerDrawGetBounds(PetscViewer,PetscInt*,const PetscReal**);
extern PetscErrorCode PetscViewerSocketSetConnection(PetscViewer,const char[],int);
extern PetscErrorCode PetscViewerBinarySkipInfo(PetscViewer);
extern PetscErrorCode PetscViewerBinarySetSkipInfo(PetscViewer,PetscBool);
extern PetscErrorCode PetscViewerBinaryGetSkipInfo(PetscViewer,PetscBool*);
extern PetscErrorCode PetscViewerBinarySetSkipOptions(PetscViewer,PetscBool );
extern PetscErrorCode PetscViewerBinaryGetSkipOptions(PetscViewer,PetscBool *);
extern PetscErrorCode PetscViewerBinarySetSkipHeader(PetscViewer,PetscBool);
extern PetscErrorCode PetscViewerBinaryGetSkipHeader(PetscViewer,PetscBool*);
extern PetscErrorCode PetscViewerBinaryReadStringArray(PetscViewer,char***);
extern PetscErrorCode PetscViewerBinaryWriteStringArray(PetscViewer,char**);

extern PetscErrorCode PetscViewerFileSetName(PetscViewer,const char[]);
extern PetscErrorCode PetscViewerFileGetName(PetscViewer,const char**);

extern PetscErrorCode PetscViewerVUGetPointer(PetscViewer, FILE**);
extern PetscErrorCode PetscViewerVUSetVecSeen(PetscViewer, PetscBool );
extern PetscErrorCode PetscViewerVUGetVecSeen(PetscViewer, PetscBool *);
extern PetscErrorCode PetscViewerVUPrintDeferred(PetscViewer, const char [], ...);
extern PetscErrorCode PetscViewerVUFlushDeferred(PetscViewer);

extern PetscErrorCode PetscViewerMathematicaInitializePackage(void);
extern PetscErrorCode PetscViewerMathematicaFinalizePackage(void);
extern PetscErrorCode PetscViewerMathematicaGetName(PetscViewer, const char **);
extern PetscErrorCode PetscViewerMathematicaSetName(PetscViewer, const char []);
extern PetscErrorCode PetscViewerMathematicaClearName(PetscViewer);
extern PetscErrorCode PetscViewerMathematicaSkipPackets(PetscViewer, int);

extern PetscErrorCode PetscViewerSiloGetName(PetscViewer, char **);
extern PetscErrorCode PetscViewerSiloSetName(PetscViewer, const char []);
extern PetscErrorCode PetscViewerSiloClearName(PetscViewer);
extern PetscErrorCode PetscViewerSiloGetMeshName(PetscViewer, char **);
extern PetscErrorCode PetscViewerSiloSetMeshName(PetscViewer, const char []);
extern PetscErrorCode PetscViewerSiloClearMeshName(PetscViewer);

extern PetscErrorCode PetscViewerNetcdfOpen(MPI_Comm,const char[],PetscFileMode,PetscViewer*);
extern PetscErrorCode PetscViewerNetcdfGetID(PetscViewer, int *);

typedef enum {PETSC_VTK_POINT_FIELD, PETSC_VTK_POINT_VECTOR_FIELD, PETSC_VTK_CELL_FIELD, PETSC_VTK_CELL_VECTOR_FIELD} PetscViewerVTKFieldType;
extern PetscErrorCode PetscViewerVTKAddField(PetscViewer,PetscObject,PetscErrorCode (*PetscViewerVTKWriteFunction)(PetscObject,PetscViewer),PetscViewerVTKFieldType,PetscObject);
extern PetscErrorCode PetscViewerVTKOpen(MPI_Comm,const char[],PetscFileMode,PetscViewer*);




extern PetscViewer PETSC_VIEWER_STDOUT_(MPI_Comm);
extern PetscErrorCode PetscViewerASCIIGetStdout(MPI_Comm,PetscViewer*);
extern PetscViewer PETSC_VIEWER_STDERR_(MPI_Comm);
extern PetscErrorCode PetscViewerASCIIGetStderr(MPI_Comm,PetscViewer*);
extern PetscViewer PETSC_VIEWER_DRAW_(MPI_Comm);
extern PetscViewer PETSC_VIEWER_SOCKET_(MPI_Comm);
extern PetscViewer PETSC_VIEWER_BINARY_(MPI_Comm);
extern PetscViewer PETSC_VIEWER_MATLAB_(MPI_Comm);
extern PetscViewer PETSC_VIEWER_HDF5_(MPI_Comm);
extern PetscViewer PETSC_VIEWER_MATHEMATICA_WORLD_PRIVATE;
# 300 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscviewer.h"
static inline PetscErrorCode PetscViewerFlowControlStart(PetscViewer viewer,PetscInt *mcnt,PetscInt *cnt)
{
  PetscErrorCode ierr;
  ;
  ierr = PetscViewerBinaryGetFlowControl(viewer,mcnt);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),304,__func__,"/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscviewer.h",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  ierr = PetscViewerBinaryGetFlowControl(viewer,cnt);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),305,__func__,"/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscviewer.h",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  return(0);
}



static inline PetscErrorCode PetscViewerFlowControlStepMaster(PetscViewer viewer,PetscInt i,PetscInt *mcnt,PetscInt cnt)
{
  PetscErrorCode ierr;
  MPI_Comm comm;

  ;
  ierr = PetscObjectGetComm((PetscObject)viewer,&comm);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),317,__func__,"/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscviewer.h",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  if (i >= *mcnt) {
    *mcnt += cnt;
    ierr = MPI_Bcast(mcnt,1,((MPI_Datatype)0x4c000405),0,comm);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),320,__func__,"/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscviewer.h",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  }
  return(0);
}



static inline PetscErrorCode PetscViewerFlowControlEndMaster(PetscViewer viewer,PetscInt *mcnt)
{
  PetscErrorCode ierr;
  MPI_Comm comm;
  ;
  ierr = PetscObjectGetComm((PetscObject)viewer,&comm);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),332,__func__,"/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscviewer.h",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  *mcnt = 0;
  ierr = MPI_Bcast(mcnt,1,((MPI_Datatype)0x4c000405),0,comm);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),334,__func__,"/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscviewer.h",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  return(0);
}



static inline PetscErrorCode PetscViewerFlowControlStepWorker(PetscViewer viewer,PetscMPIInt rank,PetscInt *mcnt)
{
  PetscErrorCode ierr;
  MPI_Comm comm;
  ;
  ierr = PetscObjectGetComm((PetscObject)viewer,&comm);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),345,__func__,"/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscviewer.h",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  while (PETSC_TRUE) {
    if (rank < *mcnt) break;
    ierr = MPI_Bcast(mcnt,1,((MPI_Datatype)0x4c000405),0,comm);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),348,__func__,"/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscviewer.h",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  }
  return(0);
}



static inline PetscErrorCode PetscViewerFlowControlEndWorker(PetscViewer viewer,PetscInt *mcnt)
{
  PetscErrorCode ierr;
  MPI_Comm comm;
  ;
  ierr = PetscObjectGetComm((PetscObject)viewer,&comm);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),360,__func__,"/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscviewer.h",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  while (PETSC_TRUE) {
    ierr = MPI_Bcast(mcnt,1,((MPI_Datatype)0x4c000405),0,comm);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),362,__func__,"/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscviewer.h",ierr,PETSC_ERROR_REPEAT," ");} while (0);
    if (!*mcnt) break;
  }
  return(0);
}




extern PetscErrorCode PetscViewerMatlabPutArray(PetscViewer,int,int,const PetscScalar*,const char*);
extern PetscErrorCode PetscViewerMatlabGetArray(PetscViewer,int,int,PetscScalar*,const char*);
extern PetscErrorCode PetscViewerMatlabPutVariable(PetscViewer,const char*,void*);
# 389 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscviewer.h"
typedef struct _n_PetscViewers* PetscViewers;
extern PetscErrorCode PetscViewersCreate(MPI_Comm,PetscViewers*);
extern PetscErrorCode PetscViewersDestroy(PetscViewers*);
extern PetscErrorCode PetscViewersGetViewer(PetscViewers,PetscInt,PetscViewer*);
# 11 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscvec.h" 2
# 21 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscvec.h"
typedef struct _p_Vec* Vec;
# 33 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscvec.h"
typedef struct _p_VecScatter* VecScatter;
# 42 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscvec.h"
typedef enum {SCATTER_FORWARD=0, SCATTER_REVERSE=1, SCATTER_FORWARD_LOCAL=2, SCATTER_REVERSE_LOCAL=3, SCATTER_LOCAL=2} ScatterMode;
# 95 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscvec.h"
typedef const char* VecType;
# 114 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscvec.h"
extern PetscClassId VEC_CLASSID;
extern PetscClassId VEC_SCATTER_CLASSID;


extern PetscErrorCode VecInitializePackage(void);
extern PetscErrorCode VecFinalizePackage(void);

extern PetscErrorCode VecCreate(MPI_Comm,Vec*);
extern PetscErrorCode VecCreateSeq(MPI_Comm,PetscInt,Vec*);
extern PetscErrorCode VecCreateMPI(MPI_Comm,PetscInt,PetscInt,Vec*);
extern PetscErrorCode VecCreateSeqWithArray(MPI_Comm,PetscInt,PetscInt,const PetscScalar[],Vec*);
extern PetscErrorCode VecCreateMPIWithArray(MPI_Comm,PetscInt,PetscInt,PetscInt,const PetscScalar[],Vec*);
extern PetscErrorCode VecCreateShared(MPI_Comm,PetscInt,PetscInt,Vec*);
extern PetscErrorCode VecSetFromOptions(Vec);
static inline PetscErrorCode VecViewFromOptions(Vec A,PetscObject B,const char name[]) {return PetscObjectViewFromOptions((PetscObject)A,B,name);}

extern PetscErrorCode VecSetUp(Vec);
extern PetscErrorCode VecDestroy(Vec*);
extern PetscErrorCode VecZeroEntries(Vec);
extern PetscErrorCode VecSetOptionsPrefix(Vec,const char[]);
extern PetscErrorCode VecAppendOptionsPrefix(Vec,const char[]);
extern PetscErrorCode VecGetOptionsPrefix(Vec,const char*[]);

extern PetscErrorCode VecSetSizes(Vec,PetscInt,PetscInt);

extern PetscErrorCode VecDotNorm2(Vec,Vec,PetscScalar*,PetscReal*);
extern PetscErrorCode VecDot(Vec,Vec,PetscScalar*);
extern PetscErrorCode VecDotRealPart(Vec,Vec,PetscReal*);
extern PetscErrorCode VecTDot(Vec,Vec,PetscScalar*);
extern PetscErrorCode VecMDot(Vec,PetscInt,const Vec[],PetscScalar[]);
extern PetscErrorCode VecMTDot(Vec,PetscInt,const Vec[],PetscScalar[]);
extern PetscErrorCode VecGetSubVector(Vec,IS,Vec*);
extern PetscErrorCode VecRestoreSubVector(Vec,IS,Vec*);
# 155 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscvec.h"
typedef enum {NORM_1=0,NORM_2=1,NORM_FROBENIUS=2,NORM_INFINITY=3,NORM_1_AND_2=4} NormType;
extern const char *const NormTypes[];
# 216 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscvec.h"
extern PetscErrorCode VecNorm(Vec,NormType,PetscReal *);
extern PetscErrorCode VecNormAvailable(Vec,NormType,PetscBool *,PetscReal *);
extern PetscErrorCode VecNormalize(Vec,PetscReal *);
extern PetscErrorCode VecSum(Vec,PetscScalar*);
extern PetscErrorCode VecMax(Vec,PetscInt*,PetscReal *);
extern PetscErrorCode VecMin(Vec,PetscInt*,PetscReal *);
extern PetscErrorCode VecScale(Vec,PetscScalar);
extern PetscErrorCode VecCopy(Vec,Vec);
extern PetscErrorCode VecSetRandom(Vec,PetscRandom);
extern PetscErrorCode VecSet(Vec,PetscScalar);
extern PetscErrorCode VecSetInf(Vec);
extern PetscErrorCode VecSwap(Vec,Vec);
extern PetscErrorCode VecAXPY(Vec,PetscScalar,Vec);
extern PetscErrorCode VecAXPBY(Vec,PetscScalar,PetscScalar,Vec);
extern PetscErrorCode VecMAXPY(Vec,PetscInt,const PetscScalar[],Vec[]);
extern PetscErrorCode VecAYPX(Vec,PetscScalar,Vec);
extern PetscErrorCode VecWAXPY(Vec,PetscScalar,Vec,Vec);
extern PetscErrorCode VecAXPBYPCZ(Vec,PetscScalar,PetscScalar,PetscScalar,Vec,Vec);
extern PetscErrorCode VecPointwiseMax(Vec,Vec,Vec);
extern PetscErrorCode VecPointwiseMaxAbs(Vec,Vec,Vec);
extern PetscErrorCode VecPointwiseMin(Vec,Vec,Vec);
extern PetscErrorCode VecPointwiseMult(Vec,Vec,Vec);
extern PetscErrorCode VecPointwiseDivide(Vec,Vec,Vec);
extern PetscErrorCode VecMaxPointwiseDivide(Vec,Vec,PetscReal*);
extern PetscErrorCode VecShift(Vec,PetscScalar);
extern PetscErrorCode VecReciprocal(Vec);
extern PetscErrorCode VecPermute(Vec, IS, PetscBool );
extern PetscErrorCode VecSqrtAbs(Vec);
extern PetscErrorCode VecLog(Vec);
extern PetscErrorCode VecExp(Vec);
extern PetscErrorCode VecAbs(Vec);
extern PetscErrorCode VecDuplicate(Vec,Vec*);
extern PetscErrorCode VecDuplicateVecs(Vec,PetscInt,Vec*[]);
extern PetscErrorCode VecDestroyVecs(PetscInt, Vec*[]);
extern PetscErrorCode VecStrideNormAll(Vec,NormType,PetscReal[]);
extern PetscErrorCode VecStrideMaxAll(Vec,PetscInt [],PetscReal []);
extern PetscErrorCode VecStrideMinAll(Vec,PetscInt [],PetscReal []);
extern PetscErrorCode VecStrideScaleAll(Vec,const PetscScalar[]);
extern PetscErrorCode VecUniqueEntries(Vec,PetscInt*,PetscScalar**);

extern PetscErrorCode VecStrideNorm(Vec,PetscInt,NormType,PetscReal*);
extern PetscErrorCode VecStrideMax(Vec,PetscInt,PetscInt *,PetscReal *);
extern PetscErrorCode VecStrideMin(Vec,PetscInt,PetscInt *,PetscReal *);
extern PetscErrorCode VecStrideScale(Vec,PetscInt,PetscScalar);
extern PetscErrorCode VecStrideSet(Vec,PetscInt,PetscScalar);


extern PetscErrorCode VecStrideGather(Vec,PetscInt,Vec,InsertMode);
extern PetscErrorCode VecStrideScatter(Vec,PetscInt,Vec,InsertMode);
extern PetscErrorCode VecStrideGatherAll(Vec,Vec[],InsertMode);
extern PetscErrorCode VecStrideScatterAll(Vec[],Vec,InsertMode);

extern PetscErrorCode VecStrideSubSetScatter(Vec,PetscInt,const PetscInt[],const PetscInt[],Vec,InsertMode);
extern PetscErrorCode VecStrideSubSetGather(Vec,PetscInt,const PetscInt[],const PetscInt[],Vec,InsertMode);

extern PetscErrorCode VecSetValues(Vec,PetscInt,const PetscInt[],const PetscScalar[],InsertMode);
extern PetscErrorCode VecGetValues(Vec,PetscInt,const PetscInt[],PetscScalar[]);
extern PetscErrorCode VecAssemblyBegin(Vec);
extern PetscErrorCode VecAssemblyEnd(Vec);
extern PetscErrorCode VecStashSetInitialSize(Vec,PetscInt,PetscInt);
extern PetscErrorCode VecStashView(Vec,PetscViewer);
extern PetscErrorCode VecStashViewFromOptions(Vec,PetscObject,const char[]);
extern PetscErrorCode VecStashGetInfo(Vec,PetscInt*,PetscInt*,PetscInt*,PetscInt*);
# 308 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscvec.h"
static inline PetscErrorCode VecSetValue(Vec v,PetscInt i,PetscScalar va,InsertMode mode) {return VecSetValues(v,1,&i,&va,mode);}


extern PetscErrorCode VecSetBlockSize(Vec,PetscInt);
extern PetscErrorCode VecGetBlockSize(Vec,PetscInt*);
extern PetscErrorCode VecSetValuesBlocked(Vec,PetscInt,const PetscInt[],const PetscScalar[],InsertMode);


extern PetscFunctionList VecList;
extern PetscErrorCode VecSetType(Vec, VecType);
extern PetscErrorCode VecGetType(Vec, VecType *);
extern PetscErrorCode VecRegister(const char[],PetscErrorCode (*)(Vec));

extern PetscErrorCode VecScatterCreate(Vec,IS,Vec,IS,VecScatter *);
extern PetscErrorCode VecScatterCreateEmpty(MPI_Comm,VecScatter *);
extern PetscErrorCode VecScatterCreateLocal(VecScatter,PetscInt,const PetscInt[],const PetscInt[],const PetscInt[],PetscInt,const PetscInt[],const PetscInt[],const PetscInt[],PetscInt);
extern PetscErrorCode VecScatterBegin(VecScatter,Vec,Vec,InsertMode,ScatterMode);
extern PetscErrorCode VecScatterEnd(VecScatter,Vec,Vec,InsertMode,ScatterMode);
extern PetscErrorCode VecScatterDestroy(VecScatter*);
extern PetscErrorCode VecScatterCopy(VecScatter,VecScatter *);
extern PetscErrorCode VecScatterView(VecScatter,PetscViewer);
static inline PetscErrorCode VecScatterViewFromOptions(VecScatter A,PetscObject obj,const char name[]) {return PetscObjectViewFromOptions((PetscObject)A,obj,name);}
extern PetscErrorCode VecScatterRemap(VecScatter,PetscInt *,PetscInt*);
extern PetscErrorCode VecScatterGetMerged(VecScatter,PetscBool *);

extern PetscErrorCode VecGetArray4d(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar****[]);
extern PetscErrorCode VecRestoreArray4d(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar****[]);
extern PetscErrorCode VecGetArray3d(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar***[]);
extern PetscErrorCode VecRestoreArray3d(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar***[]);
extern PetscErrorCode VecGetArray2d(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar**[]);
extern PetscErrorCode VecRestoreArray2d(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar**[]);
extern PetscErrorCode VecGetArray1d(Vec,PetscInt,PetscInt,PetscScalar *[]);
extern PetscErrorCode VecRestoreArray1d(Vec,PetscInt,PetscInt,PetscScalar *[]);

extern PetscErrorCode VecGetArray4dRead(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar****[]);
extern PetscErrorCode VecRestoreArray4dRead(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar****[]);
extern PetscErrorCode VecGetArray3dRead(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar***[]);
extern PetscErrorCode VecRestoreArray3dRead(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar***[]);
extern PetscErrorCode VecGetArray2dRead(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar**[]);
extern PetscErrorCode VecRestoreArray2dRead(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar**[]);
extern PetscErrorCode VecGetArray1dRead(Vec,PetscInt,PetscInt,PetscScalar *[]);
extern PetscErrorCode VecRestoreArray1dRead(Vec,PetscInt,PetscInt,PetscScalar *[]);

extern PetscErrorCode VecPlaceArray(Vec,const PetscScalar[]);
extern PetscErrorCode VecResetArray(Vec);
extern PetscErrorCode VecReplaceArray(Vec,const PetscScalar[]);
extern PetscErrorCode VecGetArrays(const Vec[],PetscInt,PetscScalar**[]);
extern PetscErrorCode VecRestoreArrays(const Vec[],PetscInt,PetscScalar**[]);

extern PetscErrorCode VecView(Vec,PetscViewer);
extern PetscErrorCode VecEqual(Vec,Vec,PetscBool *);
extern PetscErrorCode VecLoad(Vec, PetscViewer);

extern PetscErrorCode VecGetSize(Vec,PetscInt*);
extern PetscErrorCode VecGetLocalSize(Vec,PetscInt*);
extern PetscErrorCode VecGetOwnershipRange(Vec,PetscInt*,PetscInt*);
extern PetscErrorCode VecGetOwnershipRanges(Vec,const PetscInt *[]);

extern PetscErrorCode VecSetLocalToGlobalMapping(Vec,ISLocalToGlobalMapping);
extern PetscErrorCode VecSetValuesLocal(Vec,PetscInt,const PetscInt[],const PetscScalar[],InsertMode);
# 397 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscvec.h"
static inline PetscErrorCode VecSetValueLocal(Vec v,PetscInt i,PetscScalar va,InsertMode mode) {return VecSetValuesLocal(v,1,&i,&va,mode);}

extern PetscErrorCode VecSetValuesBlockedLocal(Vec,PetscInt,const PetscInt[],const PetscScalar[],InsertMode);
extern PetscErrorCode VecGetLocalToGlobalMapping(Vec,ISLocalToGlobalMapping*);

extern PetscErrorCode VecDotBegin(Vec,Vec,PetscScalar *);
extern PetscErrorCode VecDotEnd(Vec,Vec,PetscScalar *);
extern PetscErrorCode VecTDotBegin(Vec,Vec,PetscScalar *);
extern PetscErrorCode VecTDotEnd(Vec,Vec,PetscScalar *);
extern PetscErrorCode VecNormBegin(Vec,NormType,PetscReal *);
extern PetscErrorCode VecNormEnd(Vec,NormType,PetscReal *);

extern PetscErrorCode VecMDotBegin(Vec,PetscInt,const Vec[],PetscScalar[]);
extern PetscErrorCode VecMDotEnd(Vec,PetscInt,const Vec[],PetscScalar[]);
extern PetscErrorCode VecMTDotBegin(Vec,PetscInt,const Vec[],PetscScalar[]);
extern PetscErrorCode VecMTDotEnd(Vec,PetscInt,const Vec[],PetscScalar[]);
extern PetscErrorCode PetscCommSplitReductionBegin(MPI_Comm);


typedef enum {VEC_IGNORE_OFF_PROC_ENTRIES,VEC_IGNORE_NEGATIVE_INDICES} VecOption;
extern PetscErrorCode VecSetOption(Vec,VecOption,PetscBool );

extern PetscErrorCode VecGetArray(Vec,PetscScalar**);
extern PetscErrorCode VecGetArrayRead(Vec,const PetscScalar**);
extern PetscErrorCode VecRestoreArray(Vec,PetscScalar**);
extern PetscErrorCode VecRestoreArrayRead(Vec,const PetscScalar**);
extern PetscErrorCode VecGetLocalVector(Vec,Vec);
extern PetscErrorCode VecRestoreLocalVector(Vec,Vec);
extern PetscErrorCode VecGetLocalVectorRead(Vec,Vec);
extern PetscErrorCode VecRestoreLocalVectorRead(Vec,Vec);



static inline PetscErrorCode VecGetArrayPair(Vec x,Vec y,PetscScalar **xv,PetscScalar **yv)
{
  PetscErrorCode ierr;

  ;
  ierr = VecGetArray(y,yv);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),435,__func__,"/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscvec.h",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  if (x != y) {
    ierr = VecGetArrayRead(x,(const PetscScalar **)xv);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),437,__func__,"/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscvec.h",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  } else {
    *xv = *yv;
  }
  return(0);
}
static inline PetscErrorCode VecRestoreArrayPair(Vec x,Vec y,PetscScalar **xv,PetscScalar **yv)
{
  PetscErrorCode ierr;

  ;
  ierr = VecRestoreArray(y,yv);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),448,__func__,"/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscvec.h",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  if (x != y) {
    ierr = VecRestoreArrayRead(x,(const PetscScalar **)xv);do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),450,__func__,"/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscvec.h",ierr,PETSC_ERROR_REPEAT," ");} while (0);
  }
  return(0);
}
# 467 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscvec.h"
extern PetscErrorCode VecValidValues(Vec,PetscInt,PetscBool);

extern PetscErrorCode VecContourScale(Vec,PetscReal,PetscReal);





typedef enum { VECOP_VIEW = 33, VECOP_LOAD = 41, VECOP_DUPLICATE = 0} VecOperation;
extern PetscErrorCode VecSetOperation(Vec,VecOperation,void(*)(void));





extern PetscErrorCode VecMPISetGhost(Vec,PetscInt,const PetscInt[]);
extern PetscErrorCode VecCreateGhost(MPI_Comm,PetscInt,PetscInt,PetscInt,const PetscInt[],Vec*);
extern PetscErrorCode VecCreateGhostWithArray(MPI_Comm,PetscInt,PetscInt,PetscInt,const PetscInt[],const PetscScalar[],Vec*);
extern PetscErrorCode VecCreateGhostBlock(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],Vec*);
extern PetscErrorCode VecCreateGhostBlockWithArray(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],const PetscScalar[],Vec*);
extern PetscErrorCode VecGhostGetLocalForm(Vec,Vec*);
extern PetscErrorCode VecGhostRestoreLocalForm(Vec,Vec*);
extern PetscErrorCode VecGhostIsLocalForm(Vec,Vec,PetscBool*);
extern PetscErrorCode VecGhostUpdateBegin(Vec,InsertMode,ScatterMode);
extern PetscErrorCode VecGhostUpdateEnd(Vec,InsertMode,ScatterMode);

extern PetscErrorCode VecConjugate(Vec);

extern PetscErrorCode VecScatterCreateToAll(Vec,VecScatter*,Vec*);
extern PetscErrorCode VecScatterCreateToZero(Vec,VecScatter*,Vec*);

extern PetscErrorCode ISComplementVec(IS,Vec,IS*);
extern PetscErrorCode VecPow(Vec, PetscScalar);
extern PetscErrorCode VecMedian(Vec, Vec, Vec, Vec);
extern PetscErrorCode VecWhichBetween(Vec, Vec, Vec, IS *);
extern PetscErrorCode VecWhichBetweenOrEqual(Vec, Vec, Vec, IS *);
extern PetscErrorCode VecWhichGreaterThan(Vec, Vec, IS * );
extern PetscErrorCode VecWhichLessThan(Vec, Vec, IS *);
extern PetscErrorCode VecWhichEqual(Vec, Vec, IS *);
extern PetscErrorCode VecISAXPY(Vec, IS, PetscScalar,Vec);
extern PetscErrorCode VecISSet(Vec,IS, PetscScalar);
extern PetscErrorCode VecBoundGradientProjection(Vec, Vec, Vec, Vec, Vec);
extern PetscErrorCode VecStepBoundInfo(Vec,Vec,Vec,Vec,PetscReal*, PetscReal*,PetscReal*);
extern PetscErrorCode VecStepMax(Vec, Vec, PetscReal *);
extern PetscErrorCode VecStepMaxBounded(Vec,Vec,Vec,Vec,PetscReal*);

extern PetscErrorCode PetscViewerMathematicaGetVector(PetscViewer, Vec);
extern PetscErrorCode PetscViewerMathematicaPutVector(PetscViewer, Vec);
# 531 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscvec.h"
        struct _n_Vecs {PetscInt n; Vec v;};
typedef struct _n_Vecs* Vecs;
extern PetscErrorCode VecsDestroy(Vecs);
extern PetscErrorCode VecsCreateSeq(MPI_Comm,PetscInt,PetscInt,Vecs*);
extern PetscErrorCode VecsCreateSeqWithArray(MPI_Comm,PetscInt,PetscInt,PetscScalar*,Vecs*);
extern PetscErrorCode VecsDuplicate(Vecs,Vecs*);
# 562 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscvec.h"
extern PetscErrorCode VecNestGetSubVecs(Vec,PetscInt*,Vec**);
extern PetscErrorCode VecNestGetSubVec(Vec,PetscInt,Vec*);
extern PetscErrorCode VecNestSetSubVecs(Vec,PetscInt,PetscInt*,Vec*);
extern PetscErrorCode VecNestSetSubVec(Vec,PetscInt,Vec);
extern PetscErrorCode VecCreateNest(MPI_Comm,PetscInt,IS*,Vec*,Vec*);
extern PetscErrorCode VecNestGetSize(Vec,PetscInt*);

extern PetscErrorCode PetscOptionsGetVec(const char[],const char[],Vec,PetscBool*);
extern PetscErrorCode VecChop(Vec,PetscReal);

extern PetscErrorCode VecGetLayout(Vec,PetscLayout*);
extern PetscErrorCode VecSetLayout(Vec,PetscLayout);

extern PetscErrorCode PetscSectionVecView(PetscSection, Vec, PetscViewer);
extern PetscErrorCode VecGetValuesSection(Vec, PetscSection, PetscInt, PetscScalar **);
extern PetscErrorCode VecSetValuesSection(Vec, PetscSection, PetscInt, PetscScalar [], InsertMode);
extern PetscErrorCode PetscSectionVecNorm(PetscSection, PetscSection, Vec, NormType, PetscReal []);
# 7 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscpf.h" 2





extern PetscFunctionList PFList;
# 21 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscpf.h"
typedef const char* PFType;
# 38 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscpf.h"
typedef struct _p_PF* PF;

extern PetscClassId PF_CLASSID;

extern PetscErrorCode PFCreate(MPI_Comm,PetscInt,PetscInt,PF*);
extern PetscErrorCode PFSetType(PF,PFType,void*);
extern PetscErrorCode PFSet(PF,PetscErrorCode(*)(void*,PetscInt,const PetscScalar*,PetscScalar*),PetscErrorCode(*)(void*,Vec,Vec),PetscErrorCode(*)(void*,PetscViewer),PetscErrorCode(*)(void*),void*);
extern PetscErrorCode PFApply(PF,PetscInt,const PetscScalar*,PetscScalar*);
extern PetscErrorCode PFApplyVec(PF,Vec,Vec);

extern PetscErrorCode PFInitializePackage(void);

extern PetscErrorCode PFRegister(const char[],PetscErrorCode (*)(PF,void*));

extern PetscErrorCode PFDestroy(PF*);
extern PetscErrorCode PFSetFromOptions(PF);
extern PetscErrorCode PFGetType(PF,PFType*);

extern PetscErrorCode PFView(PF,PetscViewer);
static inline PetscErrorCode PFViewFromOptions(PF A,PetscObject obj,const char name[]) {return PetscObjectViewFromOptions((PetscObject)A,obj,name);}
# 4 "interpolation.c" 2
# 1 "variables.h" 1

# 1 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscdmda.h" 1



# 1 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscdm.h" 1





# 1 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmat.h" 1
# 18 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmat.h"
typedef struct _p_Mat* Mat;
# 27 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmat.h"
typedef const char* MatType;
# 129 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmat.h"
typedef enum {MAT_FACTOR_NONE, MAT_FACTOR_LU, MAT_FACTOR_CHOLESKY, MAT_FACTOR_ILU, MAT_FACTOR_ICC,MAT_FACTOR_ILUDT} MatFactorType;
extern const char *const MatFactorTypes[];

extern PetscErrorCode MatGetFactor(Mat,const char*,MatFactorType,Mat*);
extern PetscErrorCode MatGetFactorAvailable(Mat,const char*,MatFactorType,PetscBool *);
extern PetscErrorCode MatFactorGetSolverPackage(Mat,const char**);
extern PetscErrorCode MatGetFactorType(Mat,MatFactorType*);
extern PetscErrorCode MatSolverPackageRegister(const char*,const MatType,MatFactorType,PetscErrorCode(*)(Mat,MatFactorType,Mat*));
extern PetscErrorCode MatSolverPackageGet(const char*,const MatType,MatFactorType,PetscBool*,PetscBool*,PetscErrorCode (**)(Mat,MatFactorType,Mat*));



extern PetscClassId MAT_CLASSID;
extern PetscClassId MAT_COLORING_CLASSID;
extern PetscClassId MAT_FDCOLORING_CLASSID;
extern PetscClassId MAT_TRANSPOSECOLORING_CLASSID;
extern PetscClassId MAT_PARTITIONING_CLASSID;
extern PetscClassId MAT_COARSEN_CLASSID;
extern PetscClassId MAT_NULLSPACE_CLASSID;
extern PetscClassId MATMFFD_CLASSID;
# 161 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmat.h"
typedef enum {MAT_INITIAL_MATRIX,MAT_REUSE_MATRIX,MAT_IGNORE_MATRIX} MatReuse;
# 171 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmat.h"
typedef enum {MAT_DO_NOT_GET_VALUES,MAT_GET_VALUES} MatGetSubMatrixOption;

extern PetscErrorCode MatInitializePackage(void);

extern PetscErrorCode MatCreate(MPI_Comm,Mat*);
extern PetscErrorCode MatSetSizes(Mat,PetscInt,PetscInt,PetscInt,PetscInt);
extern PetscErrorCode MatSetType(Mat,MatType);
extern PetscErrorCode MatSetFromOptions(Mat);
static inline PetscErrorCode MatViewFromOptions(Mat A,PetscObject obj,const char name[]) {return PetscObjectViewFromOptions((PetscObject)A,obj,name);}
extern PetscErrorCode MatRegister(const char[],PetscErrorCode(*)(Mat));
extern PetscErrorCode MatRegisterBaseName(const char[],const char[],const char[]);
extern PetscErrorCode MatSetOptionsPrefix(Mat,const char[]);
extern PetscErrorCode MatAppendOptionsPrefix(Mat,const char[]);
extern PetscErrorCode MatGetOptionsPrefix(Mat,const char*[]);
extern PetscErrorCode MatSetErrorIfFPE(Mat,PetscBool);

extern PetscFunctionList MatList;
extern PetscFunctionList MatColoringList;
extern PetscFunctionList MatPartitioningList;
extern PetscFunctionList MatCoarsenList;
# 201 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmat.h"
typedef enum {DIFFERENT_NONZERO_PATTERN,SUBSET_NONZERO_PATTERN,SAME_NONZERO_PATTERN} MatStructure;

extern PetscErrorCode MatCreateSeqDense(MPI_Comm,PetscInt,PetscInt,PetscScalar[],Mat*);
extern PetscErrorCode MatCreateDense(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar[],Mat*);
extern PetscErrorCode MatCreateSeqAIJ(MPI_Comm,PetscInt,PetscInt,PetscInt,const PetscInt[],Mat*);
extern PetscErrorCode MatCreateAIJ(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],PetscInt,const PetscInt[],Mat*);
extern PetscErrorCode MatCreateMPIAIJWithArrays(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],const PetscInt[],const PetscScalar[],Mat *);
extern PetscErrorCode MatCreateMPIAIJWithSplitArrays(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt[],PetscInt[],PetscScalar[],PetscInt[],PetscInt[],PetscScalar[],Mat*);

extern PetscErrorCode MatCreateSeqBAIJ(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],Mat*);
extern PetscErrorCode MatCreateBAIJ(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],PetscInt,const PetscInt[],Mat*);
extern PetscErrorCode MatCreateMPIBAIJWithArrays(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],const PetscInt[],const PetscScalar[],Mat*);

extern PetscErrorCode MatCreateMPIAdj(MPI_Comm,PetscInt,PetscInt,PetscInt[],PetscInt[],PetscInt[],Mat*);
extern PetscErrorCode MatCreateSeqSBAIJ(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],Mat*);

extern PetscErrorCode MatCreateSBAIJ(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],PetscInt,const PetscInt[],Mat*);
extern PetscErrorCode MatCreateMPISBAIJWithArrays(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],const PetscInt[],const PetscScalar[],Mat *);
extern PetscErrorCode MatSeqSBAIJSetPreallocationCSR(Mat,PetscInt,const PetscInt[],const PetscInt[],const PetscScalar[]);
extern PetscErrorCode MatMPISBAIJSetPreallocationCSR(Mat,PetscInt,const PetscInt[],const PetscInt[],const PetscScalar[]);
extern PetscErrorCode MatXAIJSetPreallocation(Mat,PetscInt,const PetscInt[],const PetscInt[],const PetscInt[],const PetscInt[]);

extern PetscErrorCode MatCreateShell(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,void *,Mat*);
extern PetscErrorCode MatCreateNormal(Mat,Mat*);
extern PetscErrorCode MatCreateLRC(Mat,Mat,Mat,Mat*);
extern PetscErrorCode MatCreateIS(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,ISLocalToGlobalMapping,Mat*);
extern PetscErrorCode MatCreateSeqAIJCRL(MPI_Comm,PetscInt,PetscInt,PetscInt,const PetscInt[],Mat*);
extern PetscErrorCode MatCreateMPIAIJCRL(MPI_Comm,PetscInt,PetscInt,PetscInt,const PetscInt[],PetscInt,const PetscInt[],Mat*);

extern PetscErrorCode MatCreateSeqBSTRM(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],Mat*);
extern PetscErrorCode MatCreateMPIBSTRM(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],PetscInt,const PetscInt[],Mat*);
extern PetscErrorCode MatCreateSeqSBSTRM(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],Mat*);
extern PetscErrorCode MatCreateMPISBSTRM(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],PetscInt,const PetscInt[],Mat*);

extern PetscErrorCode MatCreateScatter(MPI_Comm,VecScatter,Mat*);
extern PetscErrorCode MatScatterSetVecScatter(Mat,VecScatter);
extern PetscErrorCode MatScatterGetVecScatter(Mat,VecScatter*);
extern PetscErrorCode MatCreateBlockMat(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt*,Mat*);
extern PetscErrorCode MatCompositeAddMat(Mat,Mat);
extern PetscErrorCode MatCompositeMerge(Mat);
extern PetscErrorCode MatCreateComposite(MPI_Comm,PetscInt,const Mat*,Mat*);
typedef enum {MAT_COMPOSITE_ADDITIVE,MAT_COMPOSITE_MULTIPLICATIVE} MatCompositeType;
extern PetscErrorCode MatCompositeSetType(Mat,MatCompositeType);

extern PetscErrorCode MatCreateFFT(MPI_Comm,PetscInt,const PetscInt[],MatType,Mat*);
extern PetscErrorCode MatCreateSeqCUFFT(MPI_Comm,PetscInt,const PetscInt[],Mat*);

extern PetscErrorCode MatCreateTranspose(Mat,Mat*);
extern PetscErrorCode MatCreateHermitianTranspose(Mat,Mat*);
extern PetscErrorCode MatCreateSubMatrix(Mat,IS,IS,Mat*);
extern PetscErrorCode MatSubMatrixUpdate(Mat,Mat,IS,IS);
extern PetscErrorCode MatCreateLocalRef(Mat,IS,IS,Mat*);

extern PetscErrorCode MatPythonSetType(Mat,const char[]);

extern PetscErrorCode MatSetUp(Mat);
extern PetscErrorCode MatDestroy(Mat*);
extern PetscErrorCode MatGetNonzeroState(Mat,PetscObjectState*);

extern PetscErrorCode MatConjugate(Mat);
extern PetscErrorCode MatRealPart(Mat);
extern PetscErrorCode MatImaginaryPart(Mat);
extern PetscErrorCode MatGetDiagonalBlock(Mat,Mat*);
extern PetscErrorCode MatGetTrace(Mat,PetscScalar*);
extern PetscErrorCode MatInvertBlockDiagonal(Mat,const PetscScalar **);


extern PetscErrorCode MatSetValues(Mat,PetscInt,const PetscInt[],PetscInt,const PetscInt[],const PetscScalar[],InsertMode);
extern PetscErrorCode MatSetValuesBlocked(Mat,PetscInt,const PetscInt[],PetscInt,const PetscInt[],const PetscScalar[],InsertMode);
extern PetscErrorCode MatSetValuesRow(Mat,PetscInt,const PetscScalar[]);
extern PetscErrorCode MatSetValuesRowLocal(Mat,PetscInt,const PetscScalar[]);
extern PetscErrorCode MatSetValuesBatch(Mat,PetscInt,PetscInt,PetscInt[],const PetscScalar[]);
extern PetscErrorCode MatSetRandom(Mat,PetscRandom);
# 287 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmat.h"
typedef struct {
  PetscInt k,j,i,c;
} MatStencil;

extern PetscErrorCode MatSetValuesStencil(Mat,PetscInt,const MatStencil[],PetscInt,const MatStencil[],const PetscScalar[],InsertMode);
extern PetscErrorCode MatSetValuesBlockedStencil(Mat,PetscInt,const MatStencil[],PetscInt,const MatStencil[],const PetscScalar[],InsertMode);
extern PetscErrorCode MatSetStencil(Mat,PetscInt,const PetscInt[],const PetscInt[],PetscInt);

extern PetscErrorCode MatSetColoring(Mat,ISColoring);
extern PetscErrorCode MatSetValuesAdifor(Mat,PetscInt,void*);
# 306 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmat.h"
typedef enum {MAT_FLUSH_ASSEMBLY=1,MAT_FINAL_ASSEMBLY=0} MatAssemblyType;
extern PetscErrorCode MatAssemblyBegin(Mat,MatAssemblyType);
extern PetscErrorCode MatAssemblyEnd(Mat,MatAssemblyType);
extern PetscErrorCode MatAssembled(Mat,PetscBool *);
# 324 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmat.h"
typedef enum {MAT_OPTION_MIN = -5,
              MAT_NEW_NONZERO_LOCATION_ERR = -4,
              MAT_UNUSED_NONZERO_LOCATION_ERR = -3,
              MAT_NEW_NONZERO_ALLOCATION_ERR = -2,
              MAT_ROW_ORIENTED = -1,
              MAT_SYMMETRIC = 1,
              MAT_STRUCTURALLY_SYMMETRIC = 2,
              MAT_NEW_DIAGONALS = 3,
              MAT_IGNORE_OFF_PROC_ENTRIES = 4,
              MAT_USE_HASH_TABLE = 5,
              MAT_KEEP_NONZERO_PATTERN = 6,
              MAT_IGNORE_ZERO_ENTRIES = 7,
              MAT_USE_INODES = 8,
              MAT_HERMITIAN = 9,
              MAT_SYMMETRY_ETERNAL = 10,
              MAT_DUMMY = 11,
              MAT_IGNORE_LOWER_TRIANGULAR = 12,
              MAT_ERROR_LOWER_TRIANGULAR = 13,
              MAT_GETROW_UPPERTRIANGULAR = 14,
              MAT_SPD = 15,
              MAT_NO_OFF_PROC_ZERO_ROWS = 16,
              MAT_NO_OFF_PROC_ENTRIES = 17,
              MAT_NEW_NONZERO_LOCATIONS = 18,
              MAT_OPTION_MAX = 19} MatOption;

extern const char *MatOptions[];
extern PetscErrorCode MatSetOption(Mat,MatOption,PetscBool);
extern PetscErrorCode MatGetOption(Mat,MatOption,PetscBool*);
extern PetscErrorCode MatGetType(Mat,MatType*);

extern PetscErrorCode MatGetValues(Mat,PetscInt,const PetscInt[],PetscInt,const PetscInt[],PetscScalar[]);
extern PetscErrorCode MatGetRow(Mat,PetscInt,PetscInt *,const PetscInt *[],const PetscScalar*[]);
extern PetscErrorCode MatRestoreRow(Mat,PetscInt,PetscInt *,const PetscInt *[],const PetscScalar*[]);
extern PetscErrorCode MatGetRowUpperTriangular(Mat);
extern PetscErrorCode MatRestoreRowUpperTriangular(Mat);
extern PetscErrorCode MatGetColumn(Mat,PetscInt,PetscInt *,const PetscInt *[],const PetscScalar*[]);
extern PetscErrorCode MatRestoreColumn(Mat,PetscInt,PetscInt *,const PetscInt *[],const PetscScalar*[]);
extern PetscErrorCode MatGetColumnVector(Mat,Vec,PetscInt);
extern PetscErrorCode MatSeqAIJGetArray(Mat,PetscScalar *[]);
extern PetscErrorCode MatSeqAIJRestoreArray(Mat,PetscScalar *[]);
extern PetscErrorCode MatSeqAIJGetMaxRowNonzeros(Mat,PetscInt*);
extern PetscErrorCode MatSeqAIJSetValuesLocalFast(Mat,PetscInt,const PetscInt[],PetscInt,const PetscInt[],const PetscScalar[],InsertMode);
extern PetscErrorCode MatDenseGetArray(Mat,PetscScalar *[]);
extern PetscErrorCode MatDenseRestoreArray(Mat,PetscScalar *[]);
extern PetscErrorCode MatGetBlockSize(Mat,PetscInt *);
extern PetscErrorCode MatSetBlockSize(Mat,PetscInt);
extern PetscErrorCode MatGetBlockSizes(Mat,PetscInt *,PetscInt *);
extern PetscErrorCode MatSetBlockSizes(Mat,PetscInt,PetscInt);
extern PetscErrorCode MatSetBlockSizesFromMats(Mat,Mat,Mat);
extern PetscErrorCode MatSetNThreads(Mat,PetscInt);
extern PetscErrorCode MatGetNThreads(Mat,PetscInt*);

extern PetscErrorCode MatMult(Mat,Vec,Vec);
extern PetscErrorCode MatMultDiagonalBlock(Mat,Vec,Vec);
extern PetscErrorCode MatMultAdd(Mat,Vec,Vec,Vec);
extern PetscErrorCode MatMultTranspose(Mat,Vec,Vec);
extern PetscErrorCode MatMultHermitianTranspose(Mat,Vec,Vec);
extern PetscErrorCode MatIsTranspose(Mat,Mat,PetscReal,PetscBool *);
extern PetscErrorCode MatIsHermitianTranspose(Mat,Mat,PetscReal,PetscBool *);
extern PetscErrorCode MatMultTransposeAdd(Mat,Vec,Vec,Vec);
extern PetscErrorCode MatMultHermitianTransposeAdd(Mat,Vec,Vec,Vec);
extern PetscErrorCode MatMultConstrained(Mat,Vec,Vec);
extern PetscErrorCode MatMultTransposeConstrained(Mat,Vec,Vec);
extern PetscErrorCode MatMatSolve(Mat,Mat,Mat);
extern PetscErrorCode MatResidual(Mat,Vec,Vec,Vec);
# 404 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmat.h"
typedef enum {MAT_DO_NOT_COPY_VALUES,MAT_COPY_VALUES,MAT_SHARE_NONZERO_PATTERN} MatDuplicateOption;

extern PetscErrorCode MatConvert(Mat,MatType,MatReuse,Mat*);
extern PetscErrorCode MatDuplicate(Mat,MatDuplicateOption,Mat*);


extern PetscErrorCode MatCopy(Mat,Mat,MatStructure);
extern PetscErrorCode MatView(Mat,PetscViewer);
extern PetscErrorCode MatIsSymmetric(Mat,PetscReal,PetscBool *);
extern PetscErrorCode MatIsStructurallySymmetric(Mat,PetscBool *);
extern PetscErrorCode MatIsHermitian(Mat,PetscReal,PetscBool *);
extern PetscErrorCode MatIsSymmetricKnown(Mat,PetscBool *,PetscBool *);
extern PetscErrorCode MatIsHermitianKnown(Mat,PetscBool *,PetscBool *);
extern PetscErrorCode MatMissingDiagonal(Mat,PetscBool *,PetscInt *);
extern PetscErrorCode MatLoad(Mat, PetscViewer);

extern PetscErrorCode MatGetRowIJ(Mat,PetscInt,PetscBool ,PetscBool ,PetscInt*,const PetscInt *[],const PetscInt *[],PetscBool *);
extern PetscErrorCode MatRestoreRowIJ(Mat,PetscInt,PetscBool ,PetscBool ,PetscInt *,const PetscInt *[],const PetscInt *[],PetscBool *);
extern PetscErrorCode MatGetColumnIJ(Mat,PetscInt,PetscBool ,PetscBool ,PetscInt*,const PetscInt *[],const PetscInt *[],PetscBool *);
extern PetscErrorCode MatRestoreColumnIJ(Mat,PetscInt,PetscBool ,PetscBool ,PetscInt *,const PetscInt *[],const PetscInt *[],PetscBool *);
# 436 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmat.h"
typedef struct {
  PetscLogDouble block_size;
  PetscLogDouble nz_allocated,nz_used,nz_unneeded;
  PetscLogDouble memory;
  PetscLogDouble assemblies;
  PetscLogDouble mallocs;
  PetscLogDouble fill_ratio_given,fill_ratio_needed;
  PetscLogDouble factor_mallocs;
} MatInfo;
# 456 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmat.h"
typedef enum {MAT_LOCAL=1,MAT_GLOBAL_MAX=2,MAT_GLOBAL_SUM=3} MatInfoType;
extern PetscErrorCode MatGetInfo(Mat,MatInfoType,MatInfo*);
extern PetscErrorCode MatGetDiagonal(Mat,Vec);
extern PetscErrorCode MatGetRowMax(Mat,Vec,PetscInt[]);
extern PetscErrorCode MatGetRowMin(Mat,Vec,PetscInt[]);
extern PetscErrorCode MatGetRowMaxAbs(Mat,Vec,PetscInt[]);
extern PetscErrorCode MatGetRowMinAbs(Mat,Vec,PetscInt[]);
extern PetscErrorCode MatGetRowSum(Mat,Vec);
extern PetscErrorCode MatTranspose(Mat,MatReuse,Mat*);
extern PetscErrorCode MatHermitianTranspose(Mat,MatReuse,Mat*);
extern PetscErrorCode MatPermute(Mat,IS,IS,Mat *);
extern PetscErrorCode MatDiagonalScale(Mat,Vec,Vec);
extern PetscErrorCode MatDiagonalSet(Mat,Vec,InsertMode);
extern PetscErrorCode MatEqual(Mat,Mat,PetscBool *);
extern PetscErrorCode MatMultEqual(Mat,Mat,PetscInt,PetscBool *);
extern PetscErrorCode MatMultAddEqual(Mat,Mat,PetscInt,PetscBool *);
extern PetscErrorCode MatMultTransposeEqual(Mat,Mat,PetscInt,PetscBool *);
extern PetscErrorCode MatMultTransposeAddEqual(Mat,Mat,PetscInt,PetscBool *);

extern PetscErrorCode MatNorm(Mat,NormType,PetscReal *);
extern PetscErrorCode MatGetColumnNorms(Mat,NormType,PetscReal *);
extern PetscErrorCode MatZeroEntries(Mat);
extern PetscErrorCode MatZeroRows(Mat,PetscInt,const PetscInt [],PetscScalar,Vec,Vec);
extern PetscErrorCode MatZeroRowsIS(Mat,IS,PetscScalar,Vec,Vec);
extern PetscErrorCode MatZeroRowsStencil(Mat,PetscInt,const MatStencil [],PetscScalar,Vec,Vec);
extern PetscErrorCode MatZeroRowsColumnsStencil(Mat,PetscInt,const MatStencil[],PetscScalar,Vec,Vec);
extern PetscErrorCode MatZeroRowsColumns(Mat,PetscInt,const PetscInt [],PetscScalar,Vec,Vec);
extern PetscErrorCode MatZeroRowsColumnsIS(Mat,IS,PetscScalar,Vec,Vec);

extern PetscErrorCode MatGetSize(Mat,PetscInt*,PetscInt*);
extern PetscErrorCode MatGetLocalSize(Mat,PetscInt*,PetscInt*);
extern PetscErrorCode MatGetOwnershipRange(Mat,PetscInt*,PetscInt*);
extern PetscErrorCode MatGetOwnershipRanges(Mat,const PetscInt**);
extern PetscErrorCode MatGetOwnershipRangeColumn(Mat,PetscInt*,PetscInt*);
extern PetscErrorCode MatGetOwnershipRangesColumn(Mat,const PetscInt**);
extern PetscErrorCode MatGetOwnershipIS(Mat,IS*,IS*);

extern PetscErrorCode MatGetSubMatrices(Mat,PetscInt,const IS[],const IS[],MatReuse,Mat *[]);
extern PetscErrorCode MatGetSubMatricesMPI(Mat,PetscInt,const IS[],const IS[],MatReuse,Mat *[]);
extern PetscErrorCode MatDestroyMatrices(PetscInt,Mat *[]);
extern PetscErrorCode MatGetSubMatrix(Mat,IS,IS,MatReuse,Mat *);
extern PetscErrorCode MatGetLocalSubMatrix(Mat,IS,IS,Mat*);
extern PetscErrorCode MatRestoreLocalSubMatrix(Mat,IS,IS,Mat*);
extern PetscErrorCode MatGetSeqNonzeroStructure(Mat,Mat*);
extern PetscErrorCode MatDestroySeqNonzeroStructure(Mat*);

extern PetscErrorCode MatCreateMPIAIJSumSeqAIJ(MPI_Comm,Mat,PetscInt,PetscInt,MatReuse,Mat*);
extern PetscErrorCode MatCreateMPIAIJSumSeqAIJSymbolic(MPI_Comm,Mat,PetscInt,PetscInt,Mat*);
extern PetscErrorCode MatCreateMPIAIJSumSeqAIJNumeric(Mat,Mat);
extern PetscErrorCode MatMPIAIJGetLocalMat(Mat,MatReuse,Mat*);
extern PetscErrorCode MatMPIAIJGetLocalMatCondensed(Mat,MatReuse,IS*,IS*,Mat*);
extern PetscErrorCode MatGetBrowsOfAcols(Mat,Mat,MatReuse,IS*,IS*,Mat*);
extern PetscErrorCode MatGetGhosts(Mat, PetscInt *,const PetscInt *[]);

extern PetscErrorCode MatIncreaseOverlap(Mat,PetscInt,IS[],PetscInt);

extern PetscErrorCode MatMatMult(Mat,Mat,MatReuse,PetscReal,Mat*);
extern PetscErrorCode MatMatMultSymbolic(Mat,Mat,PetscReal,Mat*);
extern PetscErrorCode MatMatMultNumeric(Mat,Mat,Mat);

extern PetscErrorCode MatMatMatMult(Mat,Mat,Mat,MatReuse,PetscReal,Mat*);
extern PetscErrorCode MatMatMatMultSymbolic(Mat,Mat,Mat,PetscReal,Mat*);
extern PetscErrorCode MatMatMatMultNumeric(Mat,Mat,Mat,Mat);

extern PetscErrorCode MatPtAP(Mat,Mat,MatReuse,PetscReal,Mat*);
extern PetscErrorCode MatPtAPSymbolic(Mat,Mat,PetscReal,Mat*);
extern PetscErrorCode MatPtAPNumeric(Mat,Mat,Mat);
extern PetscErrorCode MatRARt(Mat,Mat,MatReuse,PetscReal,Mat*);
extern PetscErrorCode MatRARtSymbolic(Mat,Mat,PetscReal,Mat*);
extern PetscErrorCode MatRARtNumeric(Mat,Mat,Mat);

extern PetscErrorCode MatTransposeMatMult(Mat,Mat,MatReuse,PetscReal,Mat*);
extern PetscErrorCode MatTransposetMatMultSymbolic(Mat,Mat,PetscReal,Mat*);
extern PetscErrorCode MatTransposetMatMultNumeric(Mat,Mat,Mat);
extern PetscErrorCode MatMatTransposeMult(Mat,Mat,MatReuse,PetscReal,Mat*);
extern PetscErrorCode MatMatTransposeMultSymbolic(Mat,Mat,PetscReal,Mat*);
extern PetscErrorCode MatMatTransposeMultNumeric(Mat,Mat,Mat);

extern PetscErrorCode MatAXPY(Mat,PetscScalar,Mat,MatStructure);
extern PetscErrorCode MatAYPX(Mat,PetscScalar,Mat,MatStructure);

extern PetscErrorCode MatScale(Mat,PetscScalar);
extern PetscErrorCode MatShift(Mat,PetscScalar);

extern PetscErrorCode MatSetLocalToGlobalMapping(Mat,ISLocalToGlobalMapping,ISLocalToGlobalMapping);
extern PetscErrorCode MatGetLocalToGlobalMapping(Mat,ISLocalToGlobalMapping*,ISLocalToGlobalMapping*);
extern PetscErrorCode MatGetLayouts(Mat,PetscLayout*,PetscLayout*);
extern PetscErrorCode MatZeroRowsLocal(Mat,PetscInt,const PetscInt [],PetscScalar,Vec,Vec);
extern PetscErrorCode MatZeroRowsLocalIS(Mat,IS,PetscScalar,Vec,Vec);
extern PetscErrorCode MatZeroRowsColumnsLocal(Mat,PetscInt,const PetscInt [],PetscScalar,Vec,Vec);
extern PetscErrorCode MatZeroRowsColumnsLocalIS(Mat,IS,PetscScalar,Vec,Vec);
extern PetscErrorCode MatSetValuesLocal(Mat,PetscInt,const PetscInt[],PetscInt,const PetscInt[],const PetscScalar[],InsertMode);
extern PetscErrorCode MatSetValuesBlockedLocal(Mat,PetscInt,const PetscInt[],PetscInt,const PetscInt[],const PetscScalar[],InsertMode);

extern PetscErrorCode MatStashSetInitialSize(Mat,PetscInt,PetscInt);
extern PetscErrorCode MatStashGetInfo(Mat,PetscInt*,PetscInt*,PetscInt*,PetscInt*);

extern PetscErrorCode MatInterpolate(Mat,Vec,Vec);
extern PetscErrorCode MatInterpolateAdd(Mat,Vec,Vec,Vec);
extern PetscErrorCode MatRestrict(Mat,Vec,Vec);
extern PetscErrorCode MatCreateVecs(Mat,Vec*,Vec*);
__attribute((deprecated)) static inline PetscErrorCode MatGetVecs(Mat mat,Vec *x,Vec *y) {return MatCreateVecs(mat,x,y);}
extern PetscErrorCode MatCreateRedundantMatrix(Mat,PetscInt,MPI_Comm,MatReuse,Mat*);
extern PetscErrorCode MatGetMultiProcBlock(Mat,MPI_Comm,MatReuse,Mat*);
extern PetscErrorCode MatFindZeroDiagonals(Mat,IS*);
extern PetscErrorCode MatFindOffBlockDiagonalEntries(Mat,IS*);
extern PetscErrorCode MatCreateMPIMatConcatenateSeqMat(MPI_Comm,Mat,PetscInt,MatReuse,Mat*);
# 588 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmat.h"
static inline PetscErrorCode MatSetValue(Mat v,PetscInt i,PetscInt j,PetscScalar va,InsertMode mode) {return MatSetValues(v,1,&i,1,&j,&va,mode);}

static inline PetscErrorCode MatGetValue(Mat v,PetscInt i,PetscInt j,PetscScalar *va) {return MatGetValues(v,1,&i,1,&j,va);}

static inline PetscErrorCode MatSetValueLocal(Mat v,PetscInt i,PetscInt j,PetscScalar va,InsertMode mode) {return MatSetValuesLocal(v,1,&i,1,&j,&va,mode);}
# 909 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmat.h"
extern PetscErrorCode MatShellGetContext(Mat,void *);

extern PetscErrorCode MatInodeAdjustForInodes(Mat,IS*,IS*);
extern PetscErrorCode MatInodeGetInodeSizes(Mat,PetscInt *,PetscInt *[],PetscInt *);

extern PetscErrorCode MatSeqAIJSetColumnIndices(Mat,PetscInt[]);
extern PetscErrorCode MatSeqBAIJSetColumnIndices(Mat,PetscInt[]);
extern PetscErrorCode MatCreateSeqAIJWithArrays(MPI_Comm,PetscInt,PetscInt,PetscInt[],PetscInt[],PetscScalar[],Mat*);
extern PetscErrorCode MatCreateSeqBAIJWithArrays(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt[],PetscInt[],PetscScalar[],Mat*);
extern PetscErrorCode MatCreateSeqSBAIJWithArrays(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt[],PetscInt[],PetscScalar[],Mat*);
extern PetscErrorCode MatCreateSeqAIJFromTriple(MPI_Comm,PetscInt,PetscInt,PetscInt[],PetscInt[],PetscScalar[],Mat*,PetscInt,PetscBool);



extern PetscErrorCode MatSeqBAIJSetPreallocation(Mat,PetscInt,PetscInt,const PetscInt[]);
extern PetscErrorCode MatSeqSBAIJSetPreallocation(Mat,PetscInt,PetscInt,const PetscInt[]);
extern PetscErrorCode MatSeqAIJSetPreallocation(Mat,PetscInt,const PetscInt[]);

extern PetscErrorCode MatMPIBAIJSetPreallocation(Mat,PetscInt,PetscInt,const PetscInt[],PetscInt,const PetscInt[]);
extern PetscErrorCode MatMPISBAIJSetPreallocation(Mat,PetscInt,PetscInt,const PetscInt[],PetscInt,const PetscInt[]);
extern PetscErrorCode MatMPIAIJSetPreallocation(Mat,PetscInt,const PetscInt[],PetscInt,const PetscInt[]);
extern PetscErrorCode MatSeqAIJSetPreallocationCSR(Mat,const PetscInt [],const PetscInt [],const PetscScalar []);
extern PetscErrorCode MatSeqBAIJSetPreallocationCSR(Mat,PetscInt,const PetscInt[],const PetscInt[],const PetscScalar[]);
extern PetscErrorCode MatMPIAIJSetPreallocationCSR(Mat,const PetscInt[],const PetscInt[],const PetscScalar[]);
extern PetscErrorCode MatMPIBAIJSetPreallocationCSR(Mat,PetscInt,const PetscInt[],const PetscInt[],const PetscScalar[]);
extern PetscErrorCode MatMPIAdjSetPreallocation(Mat,PetscInt[],PetscInt[],PetscInt[]);
extern PetscErrorCode MatMPIDenseSetPreallocation(Mat,PetscScalar[]);
extern PetscErrorCode MatSeqDenseSetPreallocation(Mat,PetscScalar[]);
extern PetscErrorCode MatMPIAIJGetSeqAIJ(Mat,Mat*,Mat*,const PetscInt*[]);
extern PetscErrorCode MatMPIBAIJGetSeqBAIJ(Mat,Mat*,Mat*,const PetscInt*[]);
extern PetscErrorCode MatMPIAdjCreateNonemptySubcommMat(Mat,Mat*);

extern PetscErrorCode MatISSetPreallocation(Mat,PetscInt,const PetscInt[],PetscInt,const PetscInt[]);
extern PetscErrorCode MatSeqDenseSetLDA(Mat,PetscInt);
extern PetscErrorCode MatDenseGetLocalMatrix(Mat,Mat*);

extern PetscErrorCode MatStoreValues(Mat);
extern PetscErrorCode MatRetrieveValues(Mat);

extern PetscErrorCode MatDAADSetCtx(Mat,void*);

extern PetscErrorCode MatFindNonzeroRows(Mat,IS*);
# 963 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmat.h"
typedef const char* MatOrderingType;
# 974 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmat.h"
extern PetscErrorCode MatGetOrdering(Mat,MatOrderingType,IS*,IS*);
extern PetscErrorCode MatGetOrderingList(PetscFunctionList*);
extern PetscErrorCode MatOrderingRegister(const char[],PetscErrorCode(*)(Mat,MatOrderingType,IS*,IS*));
extern PetscFunctionList MatOrderingList;

extern PetscErrorCode MatReorderForNonzeroDiagonal(Mat,PetscReal,IS,IS);
extern PetscErrorCode MatCreateLaplacian(Mat,PetscReal,PetscBool,Mat*);







typedef enum {MAT_SHIFT_NONE,MAT_SHIFT_NONZERO,MAT_SHIFT_POSITIVE_DEFINITE,MAT_SHIFT_INBLOCKS} MatFactorShiftType;
extern const char *const MatFactorShiftTypes[];
extern const char *const MatFactorShiftTypesDetail[];
# 1008 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmat.h"
typedef struct {
  PetscReal diagonal_fill;
  PetscReal usedt;
  PetscReal dt;
  PetscReal dtcol;
  PetscReal dtcount;
  PetscReal fill;
  PetscReal levels;
  PetscReal pivotinblocks;

  PetscReal zeropivot;
  PetscReal shifttype;
  PetscReal shiftamount;
} MatFactorInfo;

extern PetscErrorCode MatFactorInfoInitialize(MatFactorInfo*);
extern PetscErrorCode MatCholeskyFactor(Mat,IS,const MatFactorInfo*);
extern PetscErrorCode MatCholeskyFactorSymbolic(Mat,Mat,IS,const MatFactorInfo*);
extern PetscErrorCode MatCholeskyFactorNumeric(Mat,Mat,const MatFactorInfo*);
extern PetscErrorCode MatLUFactor(Mat,IS,IS,const MatFactorInfo*);
extern PetscErrorCode MatILUFactor(Mat,IS,IS,const MatFactorInfo*);
extern PetscErrorCode MatLUFactorSymbolic(Mat,Mat,IS,IS,const MatFactorInfo*);
extern PetscErrorCode MatILUFactorSymbolic(Mat,Mat,IS,IS,const MatFactorInfo*);
extern PetscErrorCode MatICCFactorSymbolic(Mat,Mat,IS,const MatFactorInfo*);
extern PetscErrorCode MatICCFactor(Mat,IS,const MatFactorInfo*);
extern PetscErrorCode MatLUFactorNumeric(Mat,Mat,const MatFactorInfo*);
extern PetscErrorCode MatGetInertia(Mat,PetscInt*,PetscInt*,PetscInt*);
extern PetscErrorCode MatSolve(Mat,Vec,Vec);
extern PetscErrorCode MatForwardSolve(Mat,Vec,Vec);
extern PetscErrorCode MatBackwardSolve(Mat,Vec,Vec);
extern PetscErrorCode MatSolveAdd(Mat,Vec,Vec,Vec);
extern PetscErrorCode MatSolveTranspose(Mat,Vec,Vec);
extern PetscErrorCode MatSolveTransposeAdd(Mat,Vec,Vec,Vec);
extern PetscErrorCode MatSolves(Mat,Vecs,Vecs);

extern PetscErrorCode MatSetUnfactored(Mat);
# 1058 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmat.h"
typedef enum {SOR_FORWARD_SWEEP=1,SOR_BACKWARD_SWEEP=2,SOR_SYMMETRIC_SWEEP=3,
              SOR_LOCAL_FORWARD_SWEEP=4,SOR_LOCAL_BACKWARD_SWEEP=8,
              SOR_LOCAL_SYMMETRIC_SWEEP=12,SOR_ZERO_INITIAL_GUESS=16,
              SOR_EISENSTAT=32,SOR_APPLY_UPPER=64,SOR_APPLY_LOWER=128} MatSORType;
extern PetscErrorCode MatSOR(Mat,Vec,PetscReal,MatSORType,PetscReal,PetscInt,PetscInt,Vec);
# 1077 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmat.h"
typedef struct _p_MatColoring* MatColoring;
# 1086 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmat.h"
typedef const char* MatColoringType;
# 1110 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmat.h"
typedef enum {MAT_COLORING_WEIGHT_RANDOM,MAT_COLORING_WEIGHT_LEXICAL,MAT_COLORING_WEIGHT_LF,MAT_COLORING_WEIGHT_SL} MatColoringWeightType;

extern PetscErrorCode MatColoringCreate(Mat,MatColoring*);
extern PetscErrorCode MatColoringGetDegrees(Mat,PetscInt,PetscInt*);
extern PetscErrorCode MatColoringDestroy(MatColoring*);
extern PetscErrorCode MatColoringView(MatColoring,PetscViewer);
extern PetscErrorCode MatColoringSetType(MatColoring,MatColoringType);
extern PetscErrorCode MatColoringSetFromOptions(MatColoring);
extern PetscErrorCode MatColoringSetDistance(MatColoring,PetscInt);
extern PetscErrorCode MatColoringGetDistance(MatColoring,PetscInt*);
extern PetscErrorCode MatColoringSetMaxColors(MatColoring,PetscInt);
extern PetscErrorCode MatColoringGetMaxColors(MatColoring,PetscInt*);
extern PetscErrorCode MatColoringApply(MatColoring,ISColoring*);
extern PetscErrorCode MatColoringRegister(const char[],PetscErrorCode(*)(MatColoring));
extern PetscErrorCode MatColoringPatch(Mat,PetscInt,PetscInt,ISColoringValue[],ISColoring*);
extern PetscErrorCode MatColoringSetWeightType(MatColoring,MatColoringWeightType);
extern PetscErrorCode MatColoringSetWeights(MatColoring,PetscReal*,PetscInt*);
extern PetscErrorCode MatColoringCreateWeights(MatColoring,PetscReal **,PetscInt **lperm);
# 1139 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmat.h"
typedef struct _p_MatFDColoring* MatFDColoring;

extern PetscErrorCode MatFDColoringCreate(Mat,ISColoring,MatFDColoring *);
extern PetscErrorCode MatFDColoringDestroy(MatFDColoring*);
extern PetscErrorCode MatFDColoringView(MatFDColoring,PetscViewer);
extern PetscErrorCode MatFDColoringSetFunction(MatFDColoring,PetscErrorCode (*)(void),void*);
extern PetscErrorCode MatFDColoringGetFunction(MatFDColoring,PetscErrorCode (**)(void),void**);
extern PetscErrorCode MatFDColoringSetParameters(MatFDColoring,PetscReal,PetscReal);
extern PetscErrorCode MatFDColoringSetFromOptions(MatFDColoring);
extern PetscErrorCode MatFDColoringApply(Mat,MatFDColoring,Vec,void *);
extern PetscErrorCode MatFDColoringSetF(MatFDColoring,Vec);
extern PetscErrorCode MatFDColoringGetPerturbedColumns(MatFDColoring,PetscInt*,PetscInt*[]);
extern PetscErrorCode MatFDColoringSetUp(Mat,ISColoring,MatFDColoring);
extern PetscErrorCode MatFDColoringSetBlockSize(MatFDColoring,PetscInt,PetscInt);
# 1164 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmat.h"
typedef struct _p_MatTransposeColoring* MatTransposeColoring;

extern PetscErrorCode MatTransposeColoringCreate(Mat,ISColoring,MatTransposeColoring *);
extern PetscErrorCode MatTransColoringApplySpToDen(MatTransposeColoring,Mat,Mat);
extern PetscErrorCode MatTransColoringApplyDenToSp(MatTransposeColoring,Mat,Mat);
extern PetscErrorCode MatTransposeColoringDestroy(MatTransposeColoring*);
# 1185 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmat.h"
typedef struct _p_MatPartitioning* MatPartitioning;
# 1194 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmat.h"
typedef const char* MatPartitioningType;
# 1203 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmat.h"
extern PetscErrorCode MatPartitioningCreate(MPI_Comm,MatPartitioning*);
extern PetscErrorCode MatPartitioningSetType(MatPartitioning,MatPartitioningType);
extern PetscErrorCode MatPartitioningSetNParts(MatPartitioning,PetscInt);
extern PetscErrorCode MatPartitioningSetAdjacency(MatPartitioning,Mat);
extern PetscErrorCode MatPartitioningSetVertexWeights(MatPartitioning,const PetscInt[]);
extern PetscErrorCode MatPartitioningSetPartitionWeights(MatPartitioning,const PetscReal []);
extern PetscErrorCode MatPartitioningApply(MatPartitioning,IS*);
extern PetscErrorCode MatPartitioningDestroy(MatPartitioning*);

extern PetscErrorCode MatPartitioningRegister(const char[],PetscErrorCode (*)(MatPartitioning));



extern PetscErrorCode MatPartitioningView(MatPartitioning,PetscViewer);
extern PetscErrorCode MatPartitioningSetFromOptions(MatPartitioning);
extern PetscErrorCode MatPartitioningGetType(MatPartitioning,MatPartitioningType*);

extern PetscErrorCode MatPartitioningParmetisSetCoarseSequential(MatPartitioning);
extern PetscErrorCode MatPartitioningParmetisGetEdgeCut(MatPartitioning, PetscInt *);

typedef enum { MP_CHACO_MULTILEVEL=1,MP_CHACO_SPECTRAL=2,MP_CHACO_LINEAR=4,MP_CHACO_RANDOM=5,MP_CHACO_SCATTERED=6 } MPChacoGlobalType;
extern const char *const MPChacoGlobalTypes[];
typedef enum { MP_CHACO_KERNIGHAN=1,MP_CHACO_NONE=2 } MPChacoLocalType;
extern const char *const MPChacoLocalTypes[];
typedef enum { MP_CHACO_LANCZOS=0,MP_CHACO_RQI=1 } MPChacoEigenType;
extern const char *const MPChacoEigenTypes[];

extern PetscErrorCode MatPartitioningChacoSetGlobal(MatPartitioning,MPChacoGlobalType);
extern PetscErrorCode MatPartitioningChacoGetGlobal(MatPartitioning,MPChacoGlobalType*);
extern PetscErrorCode MatPartitioningChacoSetLocal(MatPartitioning,MPChacoLocalType);
extern PetscErrorCode MatPartitioningChacoGetLocal(MatPartitioning,MPChacoLocalType*);
extern PetscErrorCode MatPartitioningChacoSetCoarseLevel(MatPartitioning,PetscReal);
extern PetscErrorCode MatPartitioningChacoSetEigenSolver(MatPartitioning,MPChacoEigenType);
extern PetscErrorCode MatPartitioningChacoGetEigenSolver(MatPartitioning,MPChacoEigenType*);
extern PetscErrorCode MatPartitioningChacoSetEigenTol(MatPartitioning,PetscReal);
extern PetscErrorCode MatPartitioningChacoGetEigenTol(MatPartitioning,PetscReal*);
extern PetscErrorCode MatPartitioningChacoSetEigenNumber(MatPartitioning,PetscInt);
extern PetscErrorCode MatPartitioningChacoGetEigenNumber(MatPartitioning,PetscInt*);
# 1250 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmat.h"
extern PetscErrorCode MatPartitioningPartySetGlobal(MatPartitioning,const char*);



extern PetscErrorCode MatPartitioningPartySetLocal(MatPartitioning,const char*);
extern PetscErrorCode MatPartitioningPartySetCoarseLevel(MatPartitioning,PetscReal);
extern PetscErrorCode MatPartitioningPartySetBipart(MatPartitioning,PetscBool);
extern PetscErrorCode MatPartitioningPartySetMatchOptimization(MatPartitioning,PetscBool);

typedef enum { MP_PTSCOTCH_QUALITY,MP_PTSCOTCH_SPEED,MP_PTSCOTCH_BALANCE,MP_PTSCOTCH_SAFETY,MP_PTSCOTCH_SCALABILITY } MPPTScotchStrategyType;
extern const char *const MPPTScotchStrategyTypes[];

extern PetscErrorCode MatPartitioningPTScotchSetImbalance(MatPartitioning,PetscReal);
extern PetscErrorCode MatPartitioningPTScotchGetImbalance(MatPartitioning,PetscReal*);
extern PetscErrorCode MatPartitioningPTScotchSetStrategy(MatPartitioning,MPPTScotchStrategyType);
extern PetscErrorCode MatPartitioningPTScotchGetStrategy(MatPartitioning,MPPTScotchStrategyType*);
# 1280 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmat.h"
typedef struct _p_MatCoarsen* MatCoarsen;
# 1289 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmat.h"
typedef const char* MatCoarsenType;




typedef struct _PetscCDIntNd{
  struct _PetscCDIntNd *next;
  PetscInt gid;
}PetscCDIntNd;


typedef struct _PetscCDArrNd{
  struct _PetscCDArrNd *next;
  struct _PetscCDIntNd *array;
}PetscCDArrNd;

typedef struct _PetscCoarsenData{
  PetscCDArrNd pool_list;
  PetscCDIntNd *new_node;
  PetscInt new_left;
  PetscInt chk_sz;
  PetscCDIntNd *extra_nodes;
  PetscCDIntNd **array;
  PetscInt size;
  Mat mat;
}PetscCoarsenData;

extern PetscErrorCode MatCoarsenCreate(MPI_Comm,MatCoarsen*);
extern PetscErrorCode MatCoarsenSetType(MatCoarsen,MatCoarsenType);
extern PetscErrorCode MatCoarsenSetAdjacency(MatCoarsen,Mat);
extern PetscErrorCode MatCoarsenSetGreedyOrdering(MatCoarsen,const IS);
extern PetscErrorCode MatCoarsenSetStrictAggs(MatCoarsen,PetscBool);
extern PetscErrorCode MatCoarsenGetData( MatCoarsen, PetscCoarsenData ** );
extern PetscErrorCode MatCoarsenApply(MatCoarsen);
extern PetscErrorCode MatCoarsenDestroy(MatCoarsen*);

extern PetscErrorCode MatCoarsenRegister(const char[],PetscErrorCode (*)(MatCoarsen));



extern PetscErrorCode MatCoarsenView(MatCoarsen,PetscViewer);
extern PetscErrorCode MatCoarsenSetFromOptions(MatCoarsen);
extern PetscErrorCode MatCoarsenGetType(MatCoarsen,MatCoarsenType*);
static inline PetscErrorCode MatCoarsenViewFromOptions(MatCoarsen A,PetscObject obj,const char name[]) {return PetscObjectViewFromOptions((PetscObject)A,obj,name);}

extern PetscErrorCode MatMeshToVertexGraph(Mat,PetscInt,Mat*);
extern PetscErrorCode MatMeshToCellGraph(Mat,PetscInt,Mat*);




typedef enum { MATOP_SET_VALUES=0,
               MATOP_GET_ROW=1,
               MATOP_RESTORE_ROW=2,
               MATOP_MULT=3,
               MATOP_MULT_ADD=4,
               MATOP_MULT_TRANSPOSE=5,
               MATOP_MULT_TRANSPOSE_ADD=6,
               MATOP_SOLVE=7,
               MATOP_SOLVE_ADD=8,
               MATOP_SOLVE_TRANSPOSE=9,
               MATOP_SOLVE_TRANSPOSE_ADD=10,
               MATOP_LUFACTOR=11,
               MATOP_CHOLESKYFACTOR=12,
               MATOP_SOR=13,
               MATOP_TRANSPOSE=14,
               MATOP_GETINFO=15,
               MATOP_EQUAL=16,
               MATOP_GET_DIAGONAL=17,
               MATOP_DIAGONAL_SCALE=18,
               MATOP_NORM=19,
               MATOP_ASSEMBLY_BEGIN=20,
               MATOP_ASSEMBLY_END=21,
               MATOP_SET_OPTION=22,
               MATOP_ZERO_ENTRIES=23,
               MATOP_ZERO_ROWS=24,
               MATOP_LUFACTOR_SYMBOLIC=25,
               MATOP_LUFACTOR_NUMERIC=26,
               MATOP_CHOLESKY_FACTOR_SYMBOLIC=27,
               MATOP_CHOLESKY_FACTOR_NUMERIC=28,
               MATOP_SETUP_PREALLOCATION=29,
               MATOP_ILUFACTOR_SYMBOLIC=30,
               MATOP_ICCFACTOR_SYMBOLIC=31,


               MATOP_DUPLICATE=34,
               MATOP_FORWARD_SOLVE=35,
               MATOP_BACKWARD_SOLVE=36,
               MATOP_ILUFACTOR=37,
               MATOP_ICCFACTOR=38,
               MATOP_AXPY=39,
               MATOP_GET_SUBMATRICES=40,
               MATOP_INCREASE_OVERLAP=41,
               MATOP_GET_VALUES=42,
               MATOP_COPY=43,
               MATOP_GET_ROW_MAX=44,
               MATOP_SCALE=45,
               MATOP_SHIFT=46,
               MATOP_DIAGONAL_SET=47,
               MATOP_ZERO_ROWS_COLUMNS=48,
               MATOP_SET_RANDOM=49,
               MATOP_GET_ROW_IJ=50,
               MATOP_RESTORE_ROW_IJ=51,
               MATOP_GET_COLUMN_IJ=52,
               MATOP_RESTORE_COLUMN_IJ=53,
               MATOP_FDCOLORING_CREATE=54,
               MATOP_COLORING_PATCH=55,
               MATOP_SET_UNFACTORED=56,
               MATOP_PERMUTE=57,
               MATOP_SET_VALUES_BLOCKED=58,
               MATOP_GET_SUBMATRIX=59,
               MATOP_DESTROY=60,
               MATOP_VIEW=61,
               MATOP_CONVERT_FROM=62,
               MATOP_MATMAT_MULT=63,
               MATOP_MATMAT_MULT_SYMBOLIC=64,
               MATOP_MATMAT_MULT_NUMERIC=65,
               MATOP_SET_LOCAL_TO_GLOBAL_MAP=66,
               MATOP_SET_VALUES_LOCAL=67,
               MATOP_ZERO_ROWS_LOCAL=68,
               MATOP_GET_ROW_MAX_ABS=69,
               MATOP_GET_ROW_MIN_ABS=70,
               MATOP_CONVERT=71,
               MATOP_SET_COLORING=72,

               MATOP_SET_VALUES_ADIFOR=74,
               MATOP_FD_COLORING_APPLY=75,
               MATOP_SET_FROM_OPTIONS=76,
               MATOP_MULT_CONSTRAINED=77,
               MATOP_MULT_TRANSPOSE_CONSTRAIN=78,
               MATOP_FIND_ZERO_DIAGONALS=79,
               MATOP_MULT_MULTIPLE=80,
               MATOP_SOLVE_MULTIPLE=81,
               MATOP_GET_INERTIA=82,
               MATOP_LOAD=83,
               MATOP_IS_SYMMETRIC=84,
               MATOP_IS_HERMITIAN=85,
               MATOP_IS_STRUCTURALLY_SYMMETRIC=86,
               MATOP_SET_VALUES_BLOCKEDLOCAL=87,
               MATOP_GET_VECS=88,
               MATOP_MAT_MULT=89,
               MATOP_MAT_MULT_SYMBOLIC=90,
               MATOP_MAT_MULT_NUMERIC=91,
               MATOP_PTAP=92,
               MATOP_PTAP_SYMBOLIC=93,
               MATOP_PTAP_NUMERIC=94,
               MATOP_MAT_TRANSPOSE_MULT=95,
               MATOP_MAT_TRANSPOSE_MULT_SYMBO=96,
               MATOP_MAT_TRANSPOSE_MULT_NUMER=97,




               MATOP_CONJUGATE=102,

               MATOP_SET_VALUES_ROW=104,
               MATOP_REAL_PART=105,
               MATOP_IMAGINARY_PART=106,
               MATOP_GET_ROW_UPPER_TRIANGULAR=107,
               MATOP_RESTORE_ROW_UPPER_TRIANG=108,
               MATOP_MAT_SOLVE=109,
               MATOP_GET_REDUNDANT_MATRIX=110,
               MATOP_GET_ROW_MIN=111,
               MATOP_GET_COLUMN_VECTOR=112,
               MATOP_MISSING_DIAGONAL=113,
               MATOP_GET_SEQ_NONZERO_STRUCTUR=114,
               MATOP_CREATE=115,
               MATOP_GET_GHOSTS=116,
               MATOP_GET_LOCAL_SUB_MATRIX=117,
               MATOP_RESTORE_LOCALSUB_MATRIX=118,
               MATOP_MULT_DIAGONAL_BLOCK=119,
               MATOP_HERMITIAN_TRANSPOSE=120,
               MATOP_MULT_HERMITIAN_TRANSPOSE=121,
               MATOP_MULT_HERMITIAN_TRANS_ADD=122,
               MATOP_GET_MULTI_PROC_BLOCK=123,
               MATOP_FIND_NONZERO_ROWS=124,
               MATOP_GET_COLUMN_NORMS=125,
               MATOP_INVERT_BLOCK_DIAGONAL=126,

               MATOP_GET_SUB_MATRICES_PARALLE=128,
               MATOP_SET_VALUES_BATCH=129,
               MATOP_TRANSPOSE_MAT_MULT=130,
               MATOP_TRANSPOSE_MAT_MULT_SYMBO=131,
               MATOP_TRANSPOSE_MAT_MULT_NUMER=132,
               MATOP_TRANSPOSE_COLORING_CREAT=133,
               MATOP_TRANS_COLORING_APPLY_SPT=134,
               MATOP_TRANS_COLORING_APPLY_DEN=135,
               MATOP_RART=136,
               MATOP_RART_SYMBOLIC=137,
               MATOP_RART_NUMERIC=138,
               MATOP_SET_BLOCK_SIZES=139,
               MATOP_AYPX=140,
               MATOP_RESIDUAL=141,
               MATOP_FDCOLORING_SETUP=142,
               MATOP_MPICONCATENATESEQ=144
             } MatOperation;
extern PetscErrorCode MatHasOperation(Mat,MatOperation,PetscBool *);
extern PetscErrorCode MatShellSetOperation(Mat,MatOperation,void(*)(void));
extern PetscErrorCode MatShellGetOperation(Mat,MatOperation,void(**)(void));
extern PetscErrorCode MatShellSetContext(Mat,void*);
# 1500 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmat.h"
extern PetscErrorCode MatMPIBAIJSetHashTableFactor(Mat,PetscReal);
extern PetscErrorCode MatISGetLocalMat(Mat,Mat*);
extern PetscErrorCode MatISSetLocalMat(Mat,Mat);
extern PetscErrorCode MatISGetMPIXAIJ(Mat,MatReuse,Mat*);
# 1518 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmat.h"
typedef struct _p_MatNullSpace* MatNullSpace;

extern PetscErrorCode MatNullSpaceCreate(MPI_Comm,PetscBool ,PetscInt,const Vec[],MatNullSpace*);
extern PetscErrorCode MatNullSpaceSetFunction(MatNullSpace,PetscErrorCode (*)(MatNullSpace,Vec,void*),void*);
extern PetscErrorCode MatNullSpaceDestroy(MatNullSpace*);
extern PetscErrorCode MatNullSpaceRemove(MatNullSpace,Vec);
extern PetscErrorCode MatGetNullSpace(Mat, MatNullSpace *);
extern PetscErrorCode MatGetTransposeNullSpace(Mat, MatNullSpace *);
extern PetscErrorCode MatSetTransposeNullSpace(Mat,MatNullSpace);
extern PetscErrorCode MatSetNullSpace(Mat,MatNullSpace);
extern PetscErrorCode MatSetNearNullSpace(Mat,MatNullSpace);
extern PetscErrorCode MatGetNearNullSpace(Mat,MatNullSpace*);
extern PetscErrorCode MatNullSpaceTest(MatNullSpace,Mat,PetscBool *);
extern PetscErrorCode MatNullSpaceView(MatNullSpace,PetscViewer);
extern PetscErrorCode MatNullSpaceGetVecs(MatNullSpace,PetscBool*,PetscInt*,const Vec**);
extern PetscErrorCode MatNullSpaceCreateRigidBody(Vec,MatNullSpace*);

extern PetscErrorCode MatReorderingSeqSBAIJ(Mat,IS);
extern PetscErrorCode MatMPISBAIJSetHashTableFactor(Mat,PetscReal);
extern PetscErrorCode MatSeqSBAIJSetColumnIndices(Mat,PetscInt *);
extern PetscErrorCode MatSeqBAIJInvertBlockDiagonal(Mat);

extern PetscErrorCode MatCreateMAIJ(Mat,PetscInt,Mat*);
extern PetscErrorCode MatMAIJRedimension(Mat,PetscInt,Mat*);
extern PetscErrorCode MatMAIJGetAIJ(Mat,Mat*);

extern PetscErrorCode MatComputeExplicitOperator(Mat,Mat*);

extern PetscErrorCode MatDiagonalScaleLocal(Mat,Vec);

extern PetscErrorCode MatCreateMFFD(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,Mat*);
extern PetscErrorCode MatMFFDSetBase(Mat,Vec,Vec);
extern PetscErrorCode MatMFFDSetFunction(Mat,PetscErrorCode(*)(void*,Vec,Vec),void*);
extern PetscErrorCode MatMFFDSetFunctioni(Mat,PetscErrorCode (*)(void*,PetscInt,Vec,PetscScalar*));
extern PetscErrorCode MatMFFDSetFunctioniBase(Mat,PetscErrorCode (*)(void*,Vec));
extern PetscErrorCode MatMFFDSetHHistory(Mat,PetscScalar[],PetscInt);
extern PetscErrorCode MatMFFDResetHHistory(Mat);
extern PetscErrorCode MatMFFDSetFunctionError(Mat,PetscReal);
extern PetscErrorCode MatMFFDSetPeriod(Mat,PetscInt);
extern PetscErrorCode MatMFFDGetH(Mat,PetscScalar *);
extern PetscErrorCode MatMFFDSetOptionsPrefix(Mat,const char[]);
extern PetscErrorCode MatMFFDCheckPositivity(void*,Vec,Vec,PetscScalar*);
extern PetscErrorCode MatMFFDSetCheckh(Mat,PetscErrorCode (*)(void*,Vec,Vec,PetscScalar*),void*);
# 1574 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmat.h"
typedef struct _p_MatMFFD* MatMFFD;
# 1583 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmat.h"
typedef const char* MatMFFDType;



extern PetscErrorCode MatMFFDSetType(Mat,MatMFFDType);
extern PetscErrorCode MatMFFDRegister(const char[],PetscErrorCode (*)(MatMFFD));

extern PetscErrorCode MatMFFDDSSetUmin(Mat,PetscReal);
extern PetscErrorCode MatMFFDWPSetComputeNormU(Mat,PetscBool );

extern PetscErrorCode PetscViewerMathematicaPutMatrix(PetscViewer, PetscInt, PetscInt, PetscReal *);
extern PetscErrorCode PetscViewerMathematicaPutCSRMatrix(PetscViewer, PetscInt, PetscInt, PetscInt *, PetscInt *, PetscReal *);
# 1741 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmat.h"
extern PetscErrorCode VecScatterPetscToFFTW(Mat,Vec,Vec);
extern PetscErrorCode VecScatterFFTWToPetsc(Mat,Vec,Vec);
extern PetscErrorCode MatCreateVecsFFTW(Mat,Vec*,Vec*,Vec*);
# 1759 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscmat.h"
extern PetscErrorCode MatCreateNest(MPI_Comm,PetscInt,const IS[],PetscInt,const IS[],const Mat[],Mat*);
extern PetscErrorCode MatNestGetSize(Mat,PetscInt*,PetscInt*);
extern PetscErrorCode MatNestGetISs(Mat,IS[],IS[]);
extern PetscErrorCode MatNestGetLocalISs(Mat,IS[],IS[]);
extern PetscErrorCode MatNestGetSubMats(Mat,PetscInt*,PetscInt*,Mat***);
extern PetscErrorCode MatNestGetSubMat(Mat,PetscInt,PetscInt,Mat*);
extern PetscErrorCode MatNestSetVecType(Mat,VecType);
extern PetscErrorCode MatNestSetSubMats(Mat,PetscInt,const IS[],PetscInt,const IS[],const Mat[]);
extern PetscErrorCode MatNestSetSubMat(Mat,PetscInt,PetscInt,Mat);

extern PetscErrorCode MatChop(Mat,PetscReal);
extern PetscErrorCode MatComputeBandwidth(Mat,PetscReal,PetscInt*);

extern PetscErrorCode MatSubdomainsCreateCoalesce(Mat,PetscInt,PetscInt*,IS**);
# 7 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscdm.h" 2
# 1 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscdmtypes.h" 1
# 15 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscdmtypes.h"
typedef struct _p_DM* DM;
# 32 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscdmtypes.h"
typedef enum {DM_BOUNDARY_NONE, DM_BOUNDARY_GHOSTED, DM_BOUNDARY_MIRROR, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_TWIST} DMBoundaryType;
# 43 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscdmtypes.h"
typedef struct _p_PetscPartitioner *PetscPartitioner;
# 8 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscdm.h" 2
# 1 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscfetypes.h" 1
# 13 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscfetypes.h"
typedef struct _p_PetscSpace *PetscSpace;
# 24 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscfetypes.h"
typedef struct _p_PetscDualSpace *PetscDualSpace;
# 35 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscfetypes.h"
typedef struct _p_PetscFE *PetscFE;
# 9 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscdm.h" 2
# 1 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscdstypes.h" 1
# 13 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscdstypes.h"
typedef struct _p_PetscDS *PetscDS;
# 10 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscdm.h" 2

extern PetscErrorCode DMInitializePackage(void);

extern PetscClassId DM_CLASSID;
# 22 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscdm.h"
typedef const char* DMType;
# 34 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscdm.h"
extern const char *const DMBoundaryTypes[];
extern PetscFunctionList DMList;
extern PetscErrorCode DMCreate(MPI_Comm,DM*);
extern PetscErrorCode DMClone(DM,DM*);
extern PetscErrorCode DMSetType(DM, DMType);
extern PetscErrorCode DMGetType(DM, DMType *);
extern PetscErrorCode DMRegister(const char[],PetscErrorCode (*)(DM));
extern PetscErrorCode DMRegisterDestroy(void);

extern PetscErrorCode DMView(DM,PetscViewer);
extern PetscErrorCode DMLoad(DM,PetscViewer);
extern PetscErrorCode DMDestroy(DM*);
extern PetscErrorCode DMCreateGlobalVector(DM,Vec*);
extern PetscErrorCode DMCreateLocalVector(DM,Vec*);
extern PetscErrorCode DMGetLocalVector(DM,Vec *);
extern PetscErrorCode DMRestoreLocalVector(DM,Vec *);
extern PetscErrorCode DMGetGlobalVector(DM,Vec *);
extern PetscErrorCode DMRestoreGlobalVector(DM,Vec *);
extern PetscErrorCode DMClearGlobalVectors(DM);
extern PetscErrorCode DMGetNamedGlobalVector(DM,const char*,Vec*);
extern PetscErrorCode DMRestoreNamedGlobalVector(DM,const char*,Vec*);
extern PetscErrorCode DMGetNamedLocalVector(DM,const char*,Vec*);
extern PetscErrorCode DMRestoreNamedLocalVector(DM,const char*,Vec*);
extern PetscErrorCode DMGetLocalToGlobalMapping(DM,ISLocalToGlobalMapping*);
extern PetscErrorCode DMCreateFieldIS(DM,PetscInt*,char***,IS**);
extern PetscErrorCode DMGetBlockSize(DM,PetscInt*);
extern PetscErrorCode DMCreateColoring(DM,ISColoringType,ISColoring*);
extern PetscErrorCode DMCreateMatrix(DM,Mat*);
extern PetscErrorCode DMSetMatrixPreallocateOnly(DM,PetscBool);
extern PetscErrorCode DMCreateInterpolation(DM,DM,Mat*,Vec*);
extern PetscErrorCode DMCreateInjection(DM,DM,Mat*);
extern PetscErrorCode DMGetWorkArray(DM,PetscInt,PetscDataType,void*);
extern PetscErrorCode DMRestoreWorkArray(DM,PetscInt,PetscDataType,void*);
extern PetscErrorCode DMRefine(DM,MPI_Comm,DM*);
extern PetscErrorCode DMCoarsen(DM,MPI_Comm,DM*);
extern PetscErrorCode DMRefineHierarchy(DM,PetscInt,DM[]);
extern PetscErrorCode DMCoarsenHierarchy(DM,PetscInt,DM[]);
extern PetscErrorCode DMCoarsenHookAdd(DM,PetscErrorCode (*)(DM,DM,void*),PetscErrorCode (*)(DM,Mat,Vec,Mat,DM,void*),void*);
extern PetscErrorCode DMRefineHookAdd(DM,PetscErrorCode (*)(DM,DM,void*),PetscErrorCode (*)(DM,Mat,DM,void*),void*);
extern PetscErrorCode DMRestrict(DM,Mat,Vec,Mat,DM);
extern PetscErrorCode DMInterpolate(DM,Mat,DM);
extern PetscErrorCode DMSetFromOptions(DM);
static inline PetscErrorCode DMViewFromOptions(DM A,PetscObject obj,const char name[]) {return PetscObjectViewFromOptions((PetscObject)A,obj,name);}

extern PetscErrorCode DMSetUp(DM);
extern PetscErrorCode DMCreateInterpolationScale(DM,DM,Mat,Vec*);
extern PetscErrorCode DMCreateAggregates(DM,DM,Mat*);
extern PetscErrorCode DMGlobalToLocalHookAdd(DM,PetscErrorCode (*)(DM,Vec,InsertMode,Vec,void*),PetscErrorCode (*)(DM,Vec,InsertMode,Vec,void*),void*);
extern PetscErrorCode DMLocalToGlobalHookAdd(DM,PetscErrorCode (*)(DM,Vec,InsertMode,Vec,void*),PetscErrorCode (*)(DM,Vec,InsertMode,Vec,void*),void*);
extern PetscErrorCode DMGlobalToLocalBegin(DM,Vec,InsertMode,Vec);
extern PetscErrorCode DMGlobalToLocalEnd(DM,Vec,InsertMode,Vec);
extern PetscErrorCode DMLocalToGlobalBegin(DM,Vec,InsertMode,Vec);
extern PetscErrorCode DMLocalToGlobalEnd(DM,Vec,InsertMode,Vec);
extern PetscErrorCode DMLocalToLocalBegin(DM,Vec,InsertMode,Vec);
extern PetscErrorCode DMLocalToLocalEnd(DM,Vec,InsertMode,Vec);
extern PetscErrorCode DMConvert(DM,DMType,DM*);


extern PetscErrorCode DMGetDimension(DM,PetscInt*);
extern PetscErrorCode DMSetDimension(DM,PetscInt);
extern PetscErrorCode DMGetDimPoints(DM,PetscInt,PetscInt*,PetscInt*);


extern PetscErrorCode DMGetCoordinateDM(DM,DM*);
extern PetscErrorCode DMSetCoordinateDM(DM,DM);
extern PetscErrorCode DMGetCoordinateDim(DM,PetscInt*);
extern PetscErrorCode DMSetCoordinateDim(DM,PetscInt);
extern PetscErrorCode DMGetCoordinateSection(DM,PetscSection*);
extern PetscErrorCode DMSetCoordinateSection(DM,PetscInt,PetscSection);
extern PetscErrorCode DMGetCoordinates(DM,Vec*);
extern PetscErrorCode DMSetCoordinates(DM,Vec);
extern PetscErrorCode DMGetCoordinatesLocal(DM,Vec*);
extern PetscErrorCode DMSetCoordinatesLocal(DM,Vec);
extern PetscErrorCode DMLocatePoints(DM,Vec,IS*);
extern PetscErrorCode DMGetPeriodicity(DM,const PetscReal**,const PetscReal**,const DMBoundaryType**);
extern PetscErrorCode DMSetPeriodicity(DM,const PetscReal[],const PetscReal[],const DMBoundaryType[]);


extern PetscErrorCode DMSubDomainHookAdd(DM,PetscErrorCode (*)(DM,DM,void*),PetscErrorCode (*)(DM,VecScatter,VecScatter,DM,void*),void*);
extern PetscErrorCode DMSubDomainRestrict(DM,VecScatter,VecScatter,DM);

extern PetscErrorCode DMSetOptionsPrefix(DM,const char []);
extern PetscErrorCode DMSetVecType(DM,VecType);
extern PetscErrorCode DMGetVecType(DM,VecType*);
extern PetscErrorCode DMSetMatType(DM,MatType);
extern PetscErrorCode DMGetMatType(DM,MatType*);
extern PetscErrorCode DMSetApplicationContext(DM,void*);
extern PetscErrorCode DMSetApplicationContextDestroy(DM,PetscErrorCode (*)(void**));
extern PetscErrorCode DMGetApplicationContext(DM,void*);
extern PetscErrorCode DMSetVariableBounds(DM,PetscErrorCode (*)(DM,Vec,Vec));
extern PetscErrorCode DMHasVariableBounds(DM,PetscBool *);
extern PetscErrorCode DMHasColoring(DM,PetscBool *);
extern PetscErrorCode DMComputeVariableBounds(DM,Vec,Vec);

extern PetscErrorCode DMCreateSubDM(DM, PetscInt, PetscInt[], IS *, DM *);
extern PetscErrorCode DMCreateFieldDecomposition(DM,PetscInt*,char***,IS**,DM**);
extern PetscErrorCode DMCreateDomainDecomposition(DM,PetscInt*,char***,IS**,IS**,DM**);
extern PetscErrorCode DMCreateDomainDecompositionScatters(DM,PetscInt,DM*,VecScatter**,VecScatter**,VecScatter**);

extern PetscErrorCode DMGetRefineLevel(DM,PetscInt*);
extern PetscErrorCode DMGetCoarsenLevel(DM,PetscInt*);
extern PetscErrorCode DMFinalizePackage(void);

extern PetscErrorCode VecGetDM(Vec, DM*);
extern PetscErrorCode VecSetDM(Vec, DM);
extern PetscErrorCode MatGetDM(Mat, DM*);
extern PetscErrorCode MatSetDM(Mat, DM);

typedef struct NLF_DAAD* NLF;




extern PetscErrorCode DMPrintCellVector(PetscInt, const char [], PetscInt, const PetscScalar []);
extern PetscErrorCode DMPrintCellMatrix(PetscInt, const char [], PetscInt, PetscInt, const PetscScalar []);
extern PetscErrorCode DMPrintLocalVec(DM, const char [], PetscReal, Vec);

extern PetscErrorCode DMGetDefaultSection(DM, PetscSection *);
extern PetscErrorCode DMSetDefaultSection(DM, PetscSection);
extern PetscErrorCode DMGetDefaultConstraints(DM, PetscSection *, Mat *);
extern PetscErrorCode DMSetDefaultConstraints(DM, PetscSection, Mat);
extern PetscErrorCode DMGetDefaultGlobalSection(DM, PetscSection *);
extern PetscErrorCode DMSetDefaultGlobalSection(DM, PetscSection);
extern PetscErrorCode DMGetDefaultSF(DM, PetscSF *);
extern PetscErrorCode DMSetDefaultSF(DM, PetscSF);
extern PetscErrorCode DMCreateDefaultSF(DM, PetscSection, PetscSection);
extern PetscErrorCode DMGetPointSF(DM, PetscSF *);
extern PetscErrorCode DMSetPointSF(DM, PetscSF);

extern PetscErrorCode DMGetOutputDM(DM, DM *);
extern PetscErrorCode DMGetOutputSequenceNumber(DM, PetscInt *, PetscReal *);
extern PetscErrorCode DMSetOutputSequenceNumber(DM, PetscInt, PetscReal);
extern PetscErrorCode DMOutputSequenceLoad(DM, PetscViewer, const char *, PetscInt, PetscReal *);

extern PetscErrorCode DMGetDS(DM, PetscDS *);
extern PetscErrorCode DMSetDS(DM, PetscDS);
extern PetscErrorCode DMGetNumFields(DM, PetscInt *);
extern PetscErrorCode DMSetNumFields(DM, PetscInt);
extern PetscErrorCode DMGetField(DM, PetscInt, PetscObject *);
extern PetscErrorCode DMSetField(DM, PetscInt, PetscObject);

typedef enum {PETSC_UNIT_LENGTH, PETSC_UNIT_MASS, PETSC_UNIT_TIME, PETSC_UNIT_CURRENT, PETSC_UNIT_TEMPERATURE, PETSC_UNIT_AMOUNT, PETSC_UNIT_LUMINOSITY, NUM_PETSC_UNITS} PetscUnit;

struct _DMInterpolationInfo {
  MPI_Comm comm;
  PetscInt dim;
  PetscInt nInput;
  PetscReal *points;
  PetscInt *cells;
  PetscInt n;
  Vec coords;
  PetscInt dof;
};
typedef struct _DMInterpolationInfo *DMInterpolationInfo;

extern PetscErrorCode DMInterpolationCreate(MPI_Comm, DMInterpolationInfo *);
extern PetscErrorCode DMInterpolationSetDim(DMInterpolationInfo, PetscInt);
extern PetscErrorCode DMInterpolationGetDim(DMInterpolationInfo, PetscInt *);
extern PetscErrorCode DMInterpolationSetDof(DMInterpolationInfo, PetscInt);
extern PetscErrorCode DMInterpolationGetDof(DMInterpolationInfo, PetscInt *);
extern PetscErrorCode DMInterpolationAddPoints(DMInterpolationInfo, PetscInt, PetscReal[]);
extern PetscErrorCode DMInterpolationSetUp(DMInterpolationInfo, DM, PetscBool);
extern PetscErrorCode DMInterpolationGetCoordinates(DMInterpolationInfo, Vec *);
extern PetscErrorCode DMInterpolationGetVector(DMInterpolationInfo, Vec *);
extern PetscErrorCode DMInterpolationRestoreVector(DMInterpolationInfo, Vec *);
extern PetscErrorCode DMInterpolationEvaluate(DMInterpolationInfo, DM, Vec, Vec);
extern PetscErrorCode DMInterpolationDestroy(DMInterpolationInfo *);
# 5 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscdmda.h" 2
# 1 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscdmdatypes.h" 1
# 14 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscdmdatypes.h"
typedef enum { DMDA_STENCIL_STAR,DMDA_STENCIL_BOX } DMDAStencilType;
# 24 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscdmdatypes.h"
typedef enum { DMDA_Q0, DMDA_Q1 } DMDAInterpolationType;
# 35 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscdmdatypes.h"
typedef enum { DMDA_ELEMENT_P1, DMDA_ELEMENT_Q1 } DMDAElementType;
# 50 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscdmdatypes.h"
typedef struct {
  PetscInt dim,dof,sw;
  PetscInt mx,my,mz;
  PetscInt xs,ys,zs;
  PetscInt xm,ym,zm;
  PetscInt gxs,gys,gzs;
  PetscInt gxm,gym,gzm;
  DMBoundaryType bx,by,bz;
  DMDAStencilType st;
  DM da;
} DMDALocalInfo;
# 6 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscdmda.h" 2

# 1 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscao.h" 1
# 19 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscao.h"
typedef struct _p_AO* AO;
# 29 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscao.h"
typedef const char* AOType;






extern PetscClassId AO_CLASSID;

extern PetscErrorCode AOInitializePackage(void);

extern PetscErrorCode AOCreate(MPI_Comm,AO*);
extern PetscErrorCode AOSetIS(AO,IS,IS);
extern PetscErrorCode AOSetFromOptions(AO);

extern PetscErrorCode AOCreateBasic(MPI_Comm,PetscInt,const PetscInt[],const PetscInt[],AO*);
extern PetscErrorCode AOCreateBasicIS(IS,IS,AO*);
extern PetscErrorCode AOCreateMemoryScalable(MPI_Comm,PetscInt,const PetscInt[],const PetscInt[],AO*);
extern PetscErrorCode AOCreateMemoryScalableIS(IS,IS,AO*);
extern PetscErrorCode AOCreateMapping(MPI_Comm,PetscInt,const PetscInt[],const PetscInt[],AO*);
extern PetscErrorCode AOCreateMappingIS(IS,IS,AO*);

extern PetscErrorCode AOView(AO,PetscViewer);
static inline PetscErrorCode AOViewFromOptions(AO A,PetscObject obj,const char name[]) {return PetscObjectViewFromOptions((PetscObject)A,obj,name);}
extern PetscErrorCode AODestroy(AO*);


extern PetscFunctionList AOList;
extern PetscErrorCode AOSetType(AO, AOType);
extern PetscErrorCode AOGetType(AO, AOType *);

extern PetscErrorCode AORegister(const char [], PetscErrorCode (*)(AO));

extern PetscErrorCode AOPetscToApplication(AO,PetscInt,PetscInt[]);
extern PetscErrorCode AOApplicationToPetsc(AO,PetscInt,PetscInt[]);
extern PetscErrorCode AOPetscToApplicationIS(AO,IS);
extern PetscErrorCode AOApplicationToPetscIS(AO,IS);

extern PetscErrorCode AOPetscToApplicationPermuteInt(AO, PetscInt, PetscInt[]);
extern PetscErrorCode AOApplicationToPetscPermuteInt(AO, PetscInt, PetscInt[]);
extern PetscErrorCode AOPetscToApplicationPermuteReal(AO, PetscInt, PetscReal[]);
extern PetscErrorCode AOApplicationToPetscPermuteReal(AO, PetscInt, PetscReal[]);

extern PetscErrorCode AOMappingHasApplicationIndex(AO, PetscInt, PetscBool *);
extern PetscErrorCode AOMappingHasPetscIndex(AO, PetscInt, PetscBool *);
# 8 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscdmda.h" 2
# 1 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscfe.h" 1






# 1 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscdt.h" 1
# 16 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscdt.h"
typedef struct _p_PetscQuadrature *PetscQuadrature;

extern PetscErrorCode PetscQuadratureCreate(MPI_Comm, PetscQuadrature *);
extern PetscErrorCode PetscQuadratureDuplicate(PetscQuadrature, PetscQuadrature *);
extern PetscErrorCode PetscQuadratureGetOrder(PetscQuadrature, PetscInt*);
extern PetscErrorCode PetscQuadratureSetOrder(PetscQuadrature, PetscInt);
extern PetscErrorCode PetscQuadratureGetData(PetscQuadrature, PetscInt*, PetscInt*, const PetscReal *[], const PetscReal *[]);
extern PetscErrorCode PetscQuadratureSetData(PetscQuadrature, PetscInt, PetscInt, const PetscReal [], const PetscReal []);
extern PetscErrorCode PetscQuadratureView(PetscQuadrature, PetscViewer);
extern PetscErrorCode PetscQuadratureDestroy(PetscQuadrature *);

extern PetscErrorCode PetscQuadratureExpandComposite(PetscQuadrature, PetscInt, const PetscReal[], const PetscReal[], PetscQuadrature *);

extern PetscErrorCode PetscDTLegendreEval(PetscInt,const PetscReal*,PetscInt,const PetscInt*,PetscReal*,PetscReal*,PetscReal*);
extern PetscErrorCode PetscDTGaussQuadrature(PetscInt,PetscReal,PetscReal,PetscReal*,PetscReal*);
extern PetscErrorCode PetscDTReconstructPoly(PetscInt,PetscInt,const PetscReal*,PetscInt,const PetscReal*,PetscReal*);
extern PetscErrorCode PetscDTGaussTensorQuadrature(PetscInt,PetscInt,PetscReal,PetscReal,PetscQuadrature*);
extern PetscErrorCode PetscDTGaussJacobiQuadrature(PetscInt,PetscInt,PetscReal,PetscReal,PetscQuadrature*);
# 8 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscfe.h" 2




typedef struct {
  PetscReal v0[3];
  PetscReal J[9];
  PetscReal invJ[9];
  PetscReal detJ;
  PetscReal n[3];
  PetscInt dim;
  PetscInt dimEmbed;
} PetscFECellGeom;

extern PetscErrorCode PetscFEInitializePackage(void);

extern PetscClassId PETSCSPACE_CLASSID;
# 33 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscfe.h"
typedef const char *PetscSpaceType;



extern PetscFunctionList PetscSpaceList;
extern PetscErrorCode PetscSpaceCreate(MPI_Comm, PetscSpace *);
extern PetscErrorCode PetscSpaceDestroy(PetscSpace *);
extern PetscErrorCode PetscSpaceSetType(PetscSpace, PetscSpaceType);
extern PetscErrorCode PetscSpaceGetType(PetscSpace, PetscSpaceType *);
extern PetscErrorCode PetscSpaceSetUp(PetscSpace);
extern PetscErrorCode PetscSpaceSetFromOptions(PetscSpace);
static inline PetscErrorCode PetscSpaceViewFromOptions(PetscSpace A,PetscObject B,const char name[]) {return PetscObjectViewFromOptions((PetscObject)A,B,name);}

extern PetscErrorCode PetscSpaceView(PetscSpace,PetscViewer);
extern PetscErrorCode PetscSpaceRegister(const char [], PetscErrorCode (*)(PetscSpace));
extern PetscErrorCode PetscSpaceRegisterDestroy(void);

extern PetscErrorCode PetscSpaceGetDimension(PetscSpace, PetscInt *);
extern PetscErrorCode PetscSpaceSetOrder(PetscSpace, PetscInt);
extern PetscErrorCode PetscSpaceGetOrder(PetscSpace, PetscInt *);
extern PetscErrorCode PetscSpaceEvaluate(PetscSpace, PetscInt, const PetscReal[], PetscReal[], PetscReal[], PetscReal[]);

extern PetscErrorCode PetscSpacePolynomialSetNumVariables(PetscSpace, PetscInt);
extern PetscErrorCode PetscSpacePolynomialGetNumVariables(PetscSpace, PetscInt *);
extern PetscErrorCode PetscSpacePolynomialSetSymmetric(PetscSpace, PetscBool);
extern PetscErrorCode PetscSpacePolynomialGetSymmetric(PetscSpace, PetscBool *);
extern PetscErrorCode PetscSpacePolynomialSetTensor(PetscSpace, PetscBool);
extern PetscErrorCode PetscSpacePolynomialGetTensor(PetscSpace, PetscBool *);

extern PetscErrorCode PetscSpaceDGSetQuadrature(PetscSpace, PetscQuadrature);
extern PetscErrorCode PetscSpaceDGGetQuadrature(PetscSpace, PetscQuadrature *);

extern PetscClassId PETSCDUALSPACE_CLASSID;
# 74 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscfe.h"
typedef const char *PetscDualSpaceType;



extern PetscFunctionList PetscDualSpaceList;
extern PetscErrorCode PetscDualSpaceCreate(MPI_Comm, PetscDualSpace *);
extern PetscErrorCode PetscDualSpaceDestroy(PetscDualSpace *);
extern PetscErrorCode PetscDualSpaceDuplicate(PetscDualSpace, PetscDualSpace *);
extern PetscErrorCode PetscDualSpaceSetType(PetscDualSpace, PetscDualSpaceType);
extern PetscErrorCode PetscDualSpaceGetType(PetscDualSpace, PetscDualSpaceType *);
extern PetscErrorCode PetscDualSpaceSetUp(PetscDualSpace);
extern PetscErrorCode PetscDualSpaceSetFromOptions(PetscDualSpace);
static inline PetscErrorCode PetscDualSpaceViewFromOptions(PetscDualSpace A,PetscObject B,const char name[]) {return PetscObjectViewFromOptions((PetscObject)A,B,name);}

extern PetscErrorCode PetscDualSpaceView(PetscDualSpace,PetscViewer);
extern PetscErrorCode PetscDualSpaceRegister(const char [], PetscErrorCode (*)(PetscDualSpace));
extern PetscErrorCode PetscDualSpaceRegisterDestroy(void);

extern PetscErrorCode PetscDualSpaceGetDimension(PetscDualSpace, PetscInt *);
extern PetscErrorCode PetscDualSpaceSetOrder(PetscDualSpace, PetscInt);
extern PetscErrorCode PetscDualSpaceGetOrder(PetscDualSpace, PetscInt *);
extern PetscErrorCode PetscDualSpaceSetDM(PetscDualSpace, DM);
extern PetscErrorCode PetscDualSpaceGetDM(PetscDualSpace, DM *);
extern PetscErrorCode PetscDualSpaceGetFunctional(PetscDualSpace, PetscInt, PetscQuadrature *);
extern PetscErrorCode PetscDualSpaceCreateReferenceCell(PetscDualSpace, PetscInt, PetscBool, DM *);

extern PetscErrorCode PetscDualSpaceApply(PetscDualSpace, PetscInt, PetscFECellGeom *, PetscInt, PetscErrorCode (*)(PetscInt, const PetscReal [], PetscInt, PetscScalar *, void *), void *, PetscScalar *);

extern PetscErrorCode PetscDualSpaceLagrangeGetContinuity(PetscDualSpace, PetscBool *);
extern PetscErrorCode PetscDualSpaceLagrangeSetContinuity(PetscDualSpace, PetscBool);

extern PetscErrorCode PetscDualSpaceGetHeightSubspace(PetscDualSpace,PetscInt,PetscDualSpace *);
extern PetscErrorCode PetscDualSpaceSimpleSetDimension(PetscDualSpace, PetscInt);
extern PetscErrorCode PetscDualSpaceSimpleSetFunctional(PetscDualSpace, PetscInt, PetscQuadrature);

extern PetscClassId PETSCFE_CLASSID;
# 120 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscfe.h"
typedef const char *PetscFEType;





extern PetscFunctionList PetscFEList;
extern PetscErrorCode PetscFECreate(MPI_Comm, PetscFE *);
extern PetscErrorCode PetscFEDestroy(PetscFE *);
extern PetscErrorCode PetscFESetType(PetscFE, PetscFEType);
extern PetscErrorCode PetscFEGetType(PetscFE, PetscFEType *);
extern PetscErrorCode PetscFESetUp(PetscFE);
extern PetscErrorCode PetscFESetFromOptions(PetscFE);
static inline PetscErrorCode PetscFEViewFromOptions(PetscFE A,PetscObject B,const char name[]) {return PetscObjectViewFromOptions((PetscObject)A,B,name);}

extern PetscErrorCode PetscFEView(PetscFE,PetscViewer);
extern PetscErrorCode PetscFERegister(const char [], PetscErrorCode (*)(PetscFE));
extern PetscErrorCode PetscFERegisterDestroy(void);
extern PetscErrorCode PetscFECreateDefault(DM, PetscInt, PetscInt, PetscBool, const char [], PetscInt, PetscFE *);

extern PetscErrorCode PetscFEGetDimension(PetscFE, PetscInt *);
extern PetscErrorCode PetscFEGetSpatialDimension(PetscFE, PetscInt *);
extern PetscErrorCode PetscFESetNumComponents(PetscFE, PetscInt);
extern PetscErrorCode PetscFEGetNumComponents(PetscFE, PetscInt *);
extern PetscErrorCode PetscFEGetTileSizes(PetscFE, PetscInt *, PetscInt *, PetscInt *, PetscInt *);
extern PetscErrorCode PetscFESetTileSizes(PetscFE, PetscInt, PetscInt, PetscInt, PetscInt);
extern PetscErrorCode PetscFESetBasisSpace(PetscFE, PetscSpace);
extern PetscErrorCode PetscFEGetBasisSpace(PetscFE, PetscSpace *);
extern PetscErrorCode PetscFESetDualSpace(PetscFE, PetscDualSpace);
extern PetscErrorCode PetscFEGetDualSpace(PetscFE, PetscDualSpace *);
extern PetscErrorCode PetscFESetQuadrature(PetscFE, PetscQuadrature);
extern PetscErrorCode PetscFEGetQuadrature(PetscFE, PetscQuadrature *);
extern PetscErrorCode PetscFEGetNumDof(PetscFE, const PetscInt **);
extern PetscErrorCode PetscFEGetDefaultTabulation(PetscFE, PetscReal **, PetscReal **, PetscReal **);
extern PetscErrorCode PetscFEGetFaceTabulation(PetscFE, PetscReal **);
extern PetscErrorCode PetscFEGetTabulation(PetscFE, PetscInt, const PetscReal[], PetscReal **, PetscReal **, PetscReal **);
extern PetscErrorCode PetscFERestoreTabulation(PetscFE, PetscInt, const PetscReal[], PetscReal **, PetscReal **, PetscReal **);
extern PetscErrorCode PetscFERefine(PetscFE, PetscFE *);

extern PetscErrorCode PetscFEIntegrate(PetscFE, PetscDS, PetscInt, PetscInt, PetscFECellGeom *, const PetscScalar[], PetscDS, const PetscScalar[], PetscReal[]);
extern PetscErrorCode PetscFEIntegrateResidual(PetscFE, PetscDS, PetscInt, PetscInt, PetscFECellGeom *, const PetscScalar[], const PetscScalar[], PetscDS, const PetscScalar[], PetscScalar[]);
extern PetscErrorCode PetscFEIntegrateBdResidual(PetscFE, PetscDS, PetscInt, PetscInt, PetscFECellGeom *, const PetscScalar[], const PetscScalar[], PetscDS, const PetscScalar[], PetscScalar[]);
extern PetscErrorCode PetscFEIntegrateJacobian(PetscFE, PetscDS, PetscInt, PetscInt, PetscInt, PetscFECellGeom *, const PetscScalar[], const PetscScalar[], PetscDS, const PetscScalar[], PetscScalar[]);
extern PetscErrorCode PetscFEIntegrateBdJacobian(PetscFE, PetscDS, PetscInt, PetscInt, PetscInt, PetscFECellGeom *, const PetscScalar[], const PetscScalar[], PetscDS, const PetscScalar[], PetscScalar[]);

extern PetscErrorCode PetscFECompositeGetMapping(PetscFE, PetscInt *, const PetscReal *[], const PetscReal *[], const PetscReal *[]);

extern PetscErrorCode PetscFEOpenCLSetRealType(PetscFE, PetscDataType);
extern PetscErrorCode PetscFEOpenCLGetRealType(PetscFE, PetscDataType *);
# 9 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscdmda.h" 2
# 31 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscdmda.h"
extern PetscErrorCode DMDASetInterpolationType(DM,DMDAInterpolationType);
extern PetscErrorCode DMDAGetInterpolationType(DM,DMDAInterpolationType*);

extern PetscErrorCode DMDASetElementType(DM,DMDAElementType);
extern PetscErrorCode DMDAGetElementType(DM,DMDAElementType*);
extern PetscErrorCode DMDAGetElements(DM,PetscInt *,PetscInt *,const PetscInt*[]);
extern PetscErrorCode DMDARestoreElements(DM,PetscInt *,PetscInt *,const PetscInt*[]);

typedef enum { DMDA_X,DMDA_Y,DMDA_Z } DMDADirection;



extern PetscErrorCode DMDACreate(MPI_Comm,DM*);
extern PetscErrorCode DMDASetSizes(DM,PetscInt,PetscInt,PetscInt);
extern PetscErrorCode DMDACreate1d(MPI_Comm,DMBoundaryType,PetscInt,PetscInt,PetscInt,const PetscInt[],DM *);
extern PetscErrorCode DMDACreate2d(MPI_Comm,DMBoundaryType,DMBoundaryType,DMDAStencilType,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],const PetscInt[],DM*);
extern PetscErrorCode DMDACreate3d(MPI_Comm,DMBoundaryType,DMBoundaryType,DMBoundaryType,DMDAStencilType,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],const PetscInt[],const PetscInt[],DM*);

extern PetscErrorCode DMDAGlobalToNaturalBegin(DM,Vec,InsertMode,Vec);
extern PetscErrorCode DMDAGlobalToNaturalEnd(DM,Vec,InsertMode,Vec);
extern PetscErrorCode DMDANaturalToGlobalBegin(DM,Vec,InsertMode,Vec);
extern PetscErrorCode DMDANaturalToGlobalEnd(DM,Vec,InsertMode,Vec);
__attribute((deprecated)) static inline PetscErrorCode DMDALocalToLocalBegin(DM dm,Vec g,InsertMode mode,Vec l) {return DMLocalToLocalBegin(dm,g,mode,l);}
__attribute((deprecated)) static inline PetscErrorCode DMDALocalToLocalEnd(DM dm,Vec g,InsertMode mode,Vec l) {return DMLocalToLocalEnd(dm,g,mode,l);}
extern PetscErrorCode DMDACreateNaturalVector(DM,Vec *);

extern PetscErrorCode DMDAGetCorners(DM,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*);
extern PetscErrorCode DMDAGetGhostCorners(DM,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*);
extern PetscErrorCode DMDAGetInfo(DM,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,DMBoundaryType*,DMBoundaryType*,DMBoundaryType*,DMDAStencilType*);
extern PetscErrorCode DMDAGetProcessorSubset(DM,DMDADirection,PetscInt,MPI_Comm*);
extern PetscErrorCode DMDAGetProcessorSubsets(DM,DMDADirection,MPI_Comm*);
extern PetscErrorCode DMDAGetRay(DM,DMDADirection,PetscInt,Vec*,VecScatter*);

extern PetscErrorCode DMDAGlobalToNaturalAllCreate(DM,VecScatter*);
extern PetscErrorCode DMDANaturalAllToGlobalCreate(DM,VecScatter*);

extern PetscErrorCode DMDAGetScatter(DM,VecScatter*,VecScatter*);
extern PetscErrorCode DMDAGetNeighbors(DM,const PetscMPIInt**);

extern PetscErrorCode DMDASetAOType(DM,AOType);
extern PetscErrorCode DMDAGetAO(DM,AO*);
extern PetscErrorCode DMDASetUniformCoordinates(DM,PetscReal,PetscReal,PetscReal,PetscReal,PetscReal,PetscReal);
extern PetscErrorCode DMDAGetCoordinateArray(DM,void*);
extern PetscErrorCode DMDARestoreCoordinateArray(DM,void*);
extern PetscErrorCode DMDAGetBoundingBox(DM,PetscReal[],PetscReal[]);
extern PetscErrorCode DMDAGetLocalBoundingBox(DM,PetscReal[],PetscReal[]);
extern PetscErrorCode DMDAGetLogicalCoordinate(DM,PetscScalar,PetscScalar,PetscScalar,PetscInt*,PetscInt*,PetscInt*,PetscScalar*,PetscScalar*,PetscScalar*);

extern PetscErrorCode DMDAMapCoordsToPeriodicDomain(DM,PetscScalar*,PetscScalar*);

extern PetscErrorCode DMDAGetReducedDMDA(DM,PetscInt,DM*);

extern PetscErrorCode DMDASetFieldName(DM,PetscInt,const char[]);
extern PetscErrorCode DMDAGetFieldName(DM,PetscInt,const char**);
extern PetscErrorCode DMDASetFieldNames(DM,const char * const *);
extern PetscErrorCode DMDAGetFieldNames(DM,const char * const **);
extern PetscErrorCode DMDASetCoordinateName(DM,PetscInt,const char[]);
extern PetscErrorCode DMDAGetCoordinateName(DM,PetscInt,const char**);

extern PetscErrorCode DMDASetBoundaryType(DM,DMBoundaryType,DMBoundaryType,DMBoundaryType);
extern PetscErrorCode DMDASetDof(DM, PetscInt);
extern PetscErrorCode DMDASetOverlap(DM,PetscInt,PetscInt,PetscInt);
extern PetscErrorCode DMDAGetOverlap(DM,PetscInt*,PetscInt*,PetscInt*);
extern PetscErrorCode DMDASetNumLocalSubDomains(DM,PetscInt);
extern PetscErrorCode DMDAGetNumLocalSubDomains(DM,PetscInt*);
extern PetscErrorCode DMDAGetOffset(DM,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*);
extern PetscErrorCode DMDASetOffset(DM,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt);
extern PetscErrorCode DMDAGetNonOverlappingRegion(DM,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*);
extern PetscErrorCode DMDASetNonOverlappingRegion(DM,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt);
extern PetscErrorCode DMDASetStencilWidth(DM, PetscInt);
extern PetscErrorCode DMDASetOwnershipRanges(DM,const PetscInt[],const PetscInt[],const PetscInt[]);
extern PetscErrorCode DMDAGetOwnershipRanges(DM,const PetscInt**,const PetscInt**,const PetscInt**);
extern PetscErrorCode DMDASetNumProcs(DM, PetscInt, PetscInt, PetscInt);
extern PetscErrorCode DMDASetStencilType(DM, DMDAStencilType);

extern PetscErrorCode DMDAVecGetArray(DM,Vec,void *);
extern PetscErrorCode DMDAVecRestoreArray(DM,Vec,void *);

extern PetscErrorCode DMDAVecGetArrayDOF(DM,Vec,void *);
extern PetscErrorCode DMDAVecRestoreArrayDOF(DM,Vec,void *);

extern PetscErrorCode DMDAVecGetArrayRead(DM,Vec,void *);
extern PetscErrorCode DMDAVecRestoreArrayRead(DM,Vec,void *);

extern PetscErrorCode DMDAVecGetArrayDOFRead(DM,Vec,void *);
extern PetscErrorCode DMDAVecRestoreArrayDOFRead(DM,Vec,void *);

extern PetscErrorCode DMDASplitComm2d(MPI_Comm,PetscInt,PetscInt,PetscInt,MPI_Comm*);

extern PetscErrorCode DMDACreatePatchIS(DM,MatStencil*,MatStencil*,IS*);
# 147 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscdmda.h"
typedef struct {PetscScalar x,y;} DMDACoor2d;
# 175 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscdmda.h"
typedef struct {PetscScalar x,y,z;} DMDACoor3d;

extern PetscErrorCode DMDAGetLocalInfo(DM,DMDALocalInfo*);

extern PetscErrorCode MatRegisterDAAD(void);
extern PetscErrorCode MatCreateDAAD(DM,Mat*);
extern PetscErrorCode MatCreateSeqUSFFT(Vec,DM,Mat*);

extern PetscErrorCode DMDASetGetMatrix(DM,PetscErrorCode (*)(DM, Mat *));
extern PetscErrorCode DMDASetBlockFills(DM,const PetscInt*,const PetscInt*);
extern PetscErrorCode DMDASetRefinementFactor(DM,PetscInt,PetscInt,PetscInt);
extern PetscErrorCode DMDAGetRefinementFactor(DM,PetscInt*,PetscInt*,PetscInt*);

extern PetscErrorCode DMDAGetArray(DM,PetscBool ,void*);
extern PetscErrorCode DMDARestoreArray(DM,PetscBool ,void*);

extern PetscErrorCode DMDACreatePF(DM,PF*);


extern PetscErrorCode DMDAGetNumCells(DM, PetscInt *, PetscInt *, PetscInt *, PetscInt *);
extern PetscErrorCode DMDAGetCellPoint(DM, PetscInt, PetscInt, PetscInt, PetscInt *);
extern PetscErrorCode DMDAGetNumVertices(DM, PetscInt *, PetscInt *, PetscInt *, PetscInt *);
extern PetscErrorCode DMDAGetNumFaces(DM, PetscInt *, PetscInt *, PetscInt *, PetscInt *, PetscInt *, PetscInt *);
extern PetscErrorCode DMDAGetHeightStratum(DM, PetscInt, PetscInt *, PetscInt *);
extern PetscErrorCode DMDAGetDepthStratum(DM, PetscInt, PetscInt *, PetscInt *);
extern PetscErrorCode DMDACreateSection(DM, const PetscInt[], const PetscInt[], const PetscInt[], PetscSection *);
extern PetscErrorCode DMDAComputeCellGeometryFEM(DM, PetscInt, PetscQuadrature, PetscReal [], PetscReal [], PetscReal [], PetscReal []);
extern PetscErrorCode DMDAGetTransitiveClosure(DM, PetscInt, PetscBool, PetscInt *, PetscInt **);
extern PetscErrorCode DMDARestoreTransitiveClosure(DM, PetscInt, PetscBool, PetscInt *, PetscInt **);
extern PetscErrorCode DMDAVecGetClosure(DM, PetscSection, Vec, PetscInt, PetscInt *, PetscScalar **);
extern PetscErrorCode DMDAVecRestoreClosure(DM, PetscSection, Vec, PetscInt, PetscInt *, PetscScalar **);
extern PetscErrorCode DMDAVecSetClosure(DM, PetscSection, Vec, PetscInt, const PetscScalar *, InsertMode);
extern PetscErrorCode DMDAGetClosure(DM, PetscSection, PetscInt, PetscInt*, const PetscInt**);
extern PetscErrorCode DMDARestoreClosure(DM, PetscSection, PetscInt, PetscInt*, const PetscInt**);
extern PetscErrorCode DMDAGetClosureScalar(DM, PetscSection, PetscInt, PetscScalar*, PetscInt*, PetscScalar**);
extern PetscErrorCode DMDARestoreClosureScalar(DM, PetscSection, PetscInt, PetscScalar*, PetscInt*, PetscScalar**);
extern PetscErrorCode DMDASetClosureScalar(DM,PetscSection,PetscInt,PetscScalar*,const PetscScalar*,InsertMode);
extern PetscErrorCode DMDAConvertToCell(DM, MatStencil, PetscInt *);
extern PetscErrorCode DMDASetVertexCoordinates(DM,PetscReal,PetscReal,PetscReal,PetscReal,PetscReal,PetscReal);
extern PetscErrorCode DMDASetPreallocationCenterDimension(DM, PetscInt);
extern PetscErrorCode DMDAGetPreallocationCenterDimension(DM, PetscInt*);

extern PetscErrorCode DMDAProjectFunction(DM, PetscErrorCode (**)(PetscInt, const PetscReal [], PetscInt, PetscScalar *, void *), void **, InsertMode, Vec);
extern PetscErrorCode DMDAProjectFunctionLocal(DM, PetscErrorCode (**)(PetscInt, const PetscReal [], PetscInt, PetscScalar *, void *), void **, InsertMode, Vec);
extern PetscErrorCode DMDAComputeL2Diff(DM, PetscErrorCode (**)(PetscInt, const PetscReal [], PetscInt, PetscScalar *, void *), void **, Vec, PetscReal *);
extern PetscErrorCode DMDAComputeL2GradientDiff(DM, PetscErrorCode (**)(PetscInt, const PetscReal [], const PetscReal [], PetscInt, PetscScalar *, void *), void **, Vec, const PetscReal[], PetscReal *);
# 3 "variables.h" 2
# 1 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscksp.h" 1





# 1 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscpc.h" 1






# 1 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscpctypes.h" 1
# 15 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscpctypes.h"
typedef struct _p_PC* PC;
# 28 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscpctypes.h"
typedef const char* PCType;
# 80 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscpctypes.h"
typedef enum { PC_SIDE_DEFAULT=-1,PC_LEFT,PC_RIGHT,PC_SYMMETRIC} PCSide;

extern const char *const *const PCSides;
# 93 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscpctypes.h"
typedef enum {
              PCRICHARDSON_CONVERGED_RTOL = 2,
              PCRICHARDSON_CONVERGED_ATOL = 3,
              PCRICHARDSON_CONVERGED_ITS = 4,
              PCRICHARDSON_DIVERGED_DTOL = -4} PCRichardsonConvergedReason;
# 106 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscpctypes.h"
typedef enum { PC_JACOBI_DIAGONAL,PC_JACOBI_ROWMAX,PC_JACOBI_ROWSUM} PCJacobiType;
extern const char *const PCJacobiTypes[];
# 128 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscpctypes.h"
typedef enum {PC_ASM_BASIC = 3,PC_ASM_RESTRICT = 1,PC_ASM_INTERPOLATE = 2,PC_ASM_NONE = 0} PCASMType;
extern const char *const PCASMTypes[];
# 159 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscpctypes.h"
typedef enum {PC_GASM_BASIC = 3,PC_GASM_RESTRICT = 1,PC_GASM_INTERPOLATE = 2,PC_GASM_NONE = 0} PCGASMType;
extern const char *const PCGASMTypes[];
# 178 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscpctypes.h"
typedef enum {PC_COMPOSITE_ADDITIVE,PC_COMPOSITE_MULTIPLICATIVE,PC_COMPOSITE_SYMMETRIC_MULTIPLICATIVE,PC_COMPOSITE_SPECIAL,PC_COMPOSITE_SCHUR} PCCompositeType;
extern const char *const PCCompositeTypes[];
# 188 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscpctypes.h"
typedef enum {PC_FIELDSPLIT_SCHUR_PRE_SELF,PC_FIELDSPLIT_SCHUR_PRE_SELFP,PC_FIELDSPLIT_SCHUR_PRE_A11,PC_FIELDSPLIT_SCHUR_PRE_USER,PC_FIELDSPLIT_SCHUR_PRE_FULL} PCFieldSplitSchurPreType;
extern const char *const PCFieldSplitSchurPreTypes[];
# 198 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscpctypes.h"
typedef enum {
  PC_FIELDSPLIT_SCHUR_FACT_DIAG,
  PC_FIELDSPLIT_SCHUR_FACT_LOWER,
  PC_FIELDSPLIT_SCHUR_FACT_UPPER,
  PC_FIELDSPLIT_SCHUR_FACT_FULL
} PCFieldSplitSchurFactType;
extern const char *const PCFieldSplitSchurFactTypes[];
# 213 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscpctypes.h"
typedef enum {PC_PARMS_GLOBAL_RAS,PC_PARMS_GLOBAL_SCHUR,PC_PARMS_GLOBAL_BJ} PCPARMSGlobalType;
extern const char *const PCPARMSGlobalTypes[];







typedef enum {PC_PARMS_LOCAL_ILU0,PC_PARMS_LOCAL_ILUK,PC_PARMS_LOCAL_ILUT,PC_PARMS_LOCAL_ARMS} PCPARMSLocalType;
extern const char *const PCPARMSLocalTypes[];
# 232 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscpctypes.h"
typedef const char *PCGAMGType;




typedef const char *PCGAMGClassicalType;
# 262 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscpctypes.h"
typedef enum { PC_MG_MULTIPLICATIVE,PC_MG_ADDITIVE,PC_MG_FULL,PC_MG_KASKADE } PCMGType;
extern const char *const PCMGTypes[];
# 278 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscpctypes.h"
typedef enum { PC_MG_CYCLE_V = 1,PC_MG_CYCLE_W = 2 } PCMGCycleType;
extern const char *const PCMGCycleTypes[];
# 288 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscpctypes.h"
typedef enum { PC_EXOTIC_FACE,PC_EXOTIC_WIREBASKET } PCExoticType;
extern const char *const PCExoticTypes[];
extern PetscErrorCode PCExoticSetType(PC,PCExoticType);
# 8 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscpc.h" 2

extern PetscErrorCode PCInitializePackage(void);





extern PetscFunctionList PCList;


extern PetscClassId PC_CLASSID;

extern PetscErrorCode PCCreate(MPI_Comm,PC*);
extern PetscErrorCode PCSetType(PC,PCType);
extern PetscErrorCode PCGetType(PC,PCType*);
extern PetscErrorCode PCSetUp(PC);
extern PetscErrorCode PCGetSetUpFailedReason(PC,PetscInt*);
extern PetscErrorCode PCSetUpOnBlocks(PC);
extern PetscErrorCode PCApply(PC,Vec,Vec);
extern PetscErrorCode PCApplySymmetricLeft(PC,Vec,Vec);
extern PetscErrorCode PCApplySymmetricRight(PC,Vec,Vec);
extern PetscErrorCode PCApplyBAorAB(PC,PCSide,Vec,Vec,Vec);
extern PetscErrorCode PCApplyTranspose(PC,Vec,Vec);
extern PetscErrorCode PCApplyTransposeExists(PC,PetscBool *);
extern PetscErrorCode PCApplyBAorABTranspose(PC,PCSide,Vec,Vec,Vec);
extern PetscErrorCode PCSetReusePreconditioner(PC,PetscBool);
extern PetscErrorCode PCGetReusePreconditioner(PC,PetscBool*);
extern PetscErrorCode PCSetErrorIfFailure(PC,PetscBool);



extern PetscErrorCode PCApplyRichardson(PC,Vec,Vec,Vec,PetscReal,PetscReal,PetscReal,PetscInt,PetscBool ,PetscInt*,PCRichardsonConvergedReason*);
extern PetscErrorCode PCApplyRichardsonExists(PC,PetscBool *);
extern PetscErrorCode PCSetInitialGuessNonzero(PC,PetscBool);
extern PetscErrorCode PCGetInitialGuessNonzero(PC,PetscBool*);
extern PetscErrorCode PCSetUseAmat(PC,PetscBool);
extern PetscErrorCode PCGetUseAmat(PC,PetscBool*);


extern PetscErrorCode PCRegister(const char[],PetscErrorCode(*)(PC));

extern PetscErrorCode PCReset(PC);
extern PetscErrorCode PCDestroy(PC*);
extern PetscErrorCode PCSetFromOptions(PC);

extern PetscErrorCode PCFactorGetMatrix(PC,Mat*);
extern PetscErrorCode PCSetModifySubMatrices(PC,PetscErrorCode(*)(PC,PetscInt,const IS[],const IS[],Mat[],void*),void*);
extern PetscErrorCode PCModifySubMatrices(PC,PetscInt,const IS[],const IS[],Mat[],void*);

extern PetscErrorCode PCSetOperators(PC,Mat,Mat);
extern PetscErrorCode PCGetOperators(PC,Mat*,Mat*);
extern PetscErrorCode PCGetOperatorsSet(PC,PetscBool *,PetscBool *);

extern PetscErrorCode PCView(PC,PetscViewer);
extern PetscErrorCode PCLoad(PC,PetscViewer);
static inline PetscErrorCode PCViewFromOptions(PC A,PetscObject obj,const char name[]) {return PetscObjectViewFromOptions((PetscObject)A,obj,name);}

extern PetscErrorCode PCSetOptionsPrefix(PC,const char[]);
extern PetscErrorCode PCAppendOptionsPrefix(PC,const char[]);
extern PetscErrorCode PCGetOptionsPrefix(PC,const char*[]);

extern PetscErrorCode PCComputeExplicitOperator(PC,Mat*);





extern PetscErrorCode PCGetDiagonalScale(PC,PetscBool *);
extern PetscErrorCode PCDiagonalScaleLeft(PC,Vec,Vec);
extern PetscErrorCode PCDiagonalScaleRight(PC,Vec,Vec);
extern PetscErrorCode PCSetDiagonalScale(PC,Vec);



extern PetscErrorCode PCJacobiSetType(PC,PCJacobiType);
extern PetscErrorCode PCJacobiGetType(PC,PCJacobiType*);
extern PetscErrorCode PCJacobiSetUseAbs(PC,PetscBool);
extern PetscErrorCode PCJacobiGetUseAbs(PC,PetscBool*);
extern PetscErrorCode PCSORSetSymmetric(PC,MatSORType);
extern PetscErrorCode PCSORGetSymmetric(PC,MatSORType*);
extern PetscErrorCode PCSORSetOmega(PC,PetscReal);
extern PetscErrorCode PCSORGetOmega(PC,PetscReal*);
extern PetscErrorCode PCSORSetIterations(PC,PetscInt,PetscInt);
extern PetscErrorCode PCSORGetIterations(PC,PetscInt*,PetscInt*);

extern PetscErrorCode PCEisenstatSetOmega(PC,PetscReal);
extern PetscErrorCode PCEisenstatGetOmega(PC,PetscReal*);
extern PetscErrorCode PCEisenstatSetNoDiagonalScaling(PC,PetscBool);
extern PetscErrorCode PCEisenstatGetNoDiagonalScaling(PC,PetscBool*);

extern PetscErrorCode PCBJacobiSetTotalBlocks(PC,PetscInt,const PetscInt[]);
extern PetscErrorCode PCBJacobiGetTotalBlocks(PC,PetscInt*,const PetscInt*[]);
extern PetscErrorCode PCBJacobiSetLocalBlocks(PC,PetscInt,const PetscInt[]);
extern PetscErrorCode PCBJacobiGetLocalBlocks(PC,PetscInt*,const PetscInt*[]);

extern PetscErrorCode PCShellSetApply(PC,PetscErrorCode (*)(PC,Vec,Vec));
extern PetscErrorCode PCShellSetApplyBA(PC,PetscErrorCode (*)(PC,PCSide,Vec,Vec,Vec));
extern PetscErrorCode PCShellSetApplyTranspose(PC,PetscErrorCode (*)(PC,Vec,Vec));
extern PetscErrorCode PCShellSetSetUp(PC,PetscErrorCode (*)(PC));
extern PetscErrorCode PCShellSetApplyRichardson(PC,PetscErrorCode (*)(PC,Vec,Vec,Vec,PetscReal,PetscReal,PetscReal,PetscInt,PetscBool ,PetscInt*,PCRichardsonConvergedReason*));
extern PetscErrorCode PCShellSetView(PC,PetscErrorCode (*)(PC,PetscViewer));
extern PetscErrorCode PCShellSetDestroy(PC,PetscErrorCode (*)(PC));
extern PetscErrorCode PCShellSetContext(PC,void*);
extern PetscErrorCode PCShellGetContext(PC,void**);
extern PetscErrorCode PCShellSetName(PC,const char[]);
extern PetscErrorCode PCShellGetName(PC,const char*[]);

extern PetscErrorCode PCFactorSetZeroPivot(PC,PetscReal);

extern PetscErrorCode PCFactorSetShiftType(PC,MatFactorShiftType);
extern PetscErrorCode PCFactorSetShiftAmount(PC,PetscReal);

extern PetscErrorCode PCFactorSetMatSolverPackage(PC,const char*);
extern PetscErrorCode PCFactorGetMatSolverPackage(PC,const char**);
extern PetscErrorCode PCFactorSetUpMatSolverPackage(PC);

extern PetscErrorCode PCFactorSetFill(PC,PetscReal);
extern PetscErrorCode PCFactorSetColumnPivot(PC,PetscReal);
extern PetscErrorCode PCFactorReorderForNonzeroDiagonal(PC,PetscReal);

extern PetscErrorCode PCFactorSetMatOrderingType(PC,MatOrderingType);
extern PetscErrorCode PCFactorSetReuseOrdering(PC,PetscBool );
extern PetscErrorCode PCFactorSetReuseFill(PC,PetscBool );
extern PetscErrorCode PCFactorSetUseInPlace(PC,PetscBool);
extern PetscErrorCode PCFactorGetUseInPlace(PC,PetscBool*);
extern PetscErrorCode PCFactorSetAllowDiagonalFill(PC,PetscBool);
extern PetscErrorCode PCFactorGetAllowDiagonalFill(PC,PetscBool*);
extern PetscErrorCode PCFactorSetPivotInBlocks(PC,PetscBool);

extern PetscErrorCode PCFactorSetLevels(PC,PetscInt);
extern PetscErrorCode PCFactorGetLevels(PC,PetscInt*);
extern PetscErrorCode PCFactorSetDropTolerance(PC,PetscReal,PetscReal,PetscInt);

extern PetscErrorCode PCASMSetLocalSubdomains(PC,PetscInt,IS[],IS[]);
extern PetscErrorCode PCASMSetTotalSubdomains(PC,PetscInt,IS[],IS[]);
extern PetscErrorCode PCASMSetOverlap(PC,PetscInt);
extern PetscErrorCode PCASMSetDMSubdomains(PC,PetscBool);
extern PetscErrorCode PCASMGetDMSubdomains(PC,PetscBool*);
extern PetscErrorCode PCASMSetSortIndices(PC,PetscBool);

extern PetscErrorCode PCASMSetType(PC,PCASMType);
extern PetscErrorCode PCASMGetType(PC,PCASMType*);
extern PetscErrorCode PCASMSetLocalType(PC,PCCompositeType);
extern PetscErrorCode PCASMGetLocalType(PC,PCCompositeType*);
extern PetscErrorCode PCASMCreateSubdomains(Mat,PetscInt,IS*[]);
extern PetscErrorCode PCASMDestroySubdomains(PetscInt,IS[],IS[]);
extern PetscErrorCode PCASMCreateSubdomains2D(PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt*,IS**,IS**);
extern PetscErrorCode PCASMGetLocalSubdomains(PC,PetscInt*,IS*[],IS*[]);
extern PetscErrorCode PCASMGetLocalSubmatrices(PC,PetscInt*,Mat*[]);

extern PetscErrorCode PCGASMSetTotalSubdomains(PC,PetscInt);
extern PetscErrorCode PCGASMSetSubdomains(PC,PetscInt,IS[],IS[]);
extern PetscErrorCode PCGASMSetOverlap(PC,PetscInt);
extern PetscErrorCode PCGASMSetUseDMSubdomains(PC,PetscBool);
extern PetscErrorCode PCGASMGetUseDMSubdomains(PC,PetscBool*);
extern PetscErrorCode PCGASMSetSortIndices(PC,PetscBool );

extern PetscErrorCode PCGASMSetType(PC,PCGASMType);
extern PetscErrorCode PCGASMCreateSubdomains(Mat,PetscInt,PetscInt*,IS*[]);
extern PetscErrorCode PCGASMDestroySubdomains(PetscInt,IS*[],IS*[]);
extern PetscErrorCode PCGASMCreateSubdomains2D(PC,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt*,IS**,IS**);
extern PetscErrorCode PCGASMGetSubdomains(PC,PetscInt*,IS*[],IS*[]);
extern PetscErrorCode PCGASMGetSubmatrices(PC,PetscInt*,Mat*[]);

extern PetscErrorCode PCCompositeSetType(PC,PCCompositeType);
extern PetscErrorCode PCCompositeGetType(PC,PCCompositeType*);
extern PetscErrorCode PCCompositeAddPC(PC,PCType);
extern PetscErrorCode PCCompositeGetNumberPC(PC,PetscInt *);
extern PetscErrorCode PCCompositeGetPC(PC,PetscInt,PC *);
extern PetscErrorCode PCCompositeSpecialSetAlpha(PC,PetscScalar);

extern PetscErrorCode PCRedundantSetNumber(PC,PetscInt);
extern PetscErrorCode PCRedundantSetScatter(PC,VecScatter,VecScatter);
extern PetscErrorCode PCRedundantGetOperators(PC,Mat*,Mat*);

extern PetscErrorCode PCSPAISetEpsilon(PC,double);
extern PetscErrorCode PCSPAISetNBSteps(PC,PetscInt);
extern PetscErrorCode PCSPAISetMax(PC,PetscInt);
extern PetscErrorCode PCSPAISetMaxNew(PC,PetscInt);
extern PetscErrorCode PCSPAISetBlockSize(PC,PetscInt);
extern PetscErrorCode PCSPAISetCacheSize(PC,PetscInt);
extern PetscErrorCode PCSPAISetVerbose(PC,PetscInt);
extern PetscErrorCode PCSPAISetSp(PC,PetscInt);

extern PetscErrorCode PCHYPRESetType(PC,const char[]);
extern PetscErrorCode PCHYPREGetType(PC,const char*[]);
extern PetscErrorCode PCHYPRESetDiscreteGradient(PC,Mat);
extern PetscErrorCode PCHYPRESetDiscreteCurl(PC,Mat);
extern PetscErrorCode PCHYPRESetEdgeConstantVectors(PC,Vec,Vec,Vec);
extern PetscErrorCode PCHYPRESetAlphaPoissonMatrix(PC,Mat);
extern PetscErrorCode PCHYPRESetBetaPoissonMatrix(PC,Mat);
extern PetscErrorCode PCBJacobiGetLocalBlocks(PC,PetscInt*,const PetscInt*[]);
extern PetscErrorCode PCBJacobiGetTotalBlocks(PC,PetscInt*,const PetscInt*[]);

extern PetscErrorCode PCFieldSplitSetFields(PC,const char[],PetscInt,const PetscInt*,const PetscInt*);
extern PetscErrorCode PCFieldSplitSetType(PC,PCCompositeType);
extern PetscErrorCode PCFieldSplitGetType(PC,PCCompositeType*);
extern PetscErrorCode PCFieldSplitSetBlockSize(PC,PetscInt);
extern PetscErrorCode PCFieldSplitSetIS(PC,const char[],IS);
extern PetscErrorCode PCFieldSplitGetIS(PC,const char[],IS*);
extern PetscErrorCode PCFieldSplitSetDMSplits(PC,PetscBool);
extern PetscErrorCode PCFieldSplitGetDMSplits(PC,PetscBool*);
extern PetscErrorCode PCFieldSplitSetDiagUseAmat(PC,PetscBool);
extern PetscErrorCode PCFieldSplitGetDiagUseAmat(PC,PetscBool*);
extern PetscErrorCode PCFieldSplitSetOffDiagUseAmat(PC,PetscBool);
extern PetscErrorCode PCFieldSplitGetOffDiagUseAmat(PC,PetscBool*);

extern __attribute((deprecated)) PetscErrorCode PCFieldSplitSchurPrecondition(PC,PCFieldSplitSchurPreType,Mat);
extern PetscErrorCode PCFieldSplitSetSchurPre(PC,PCFieldSplitSchurPreType,Mat);
extern PetscErrorCode PCFieldSplitGetSchurPre(PC,PCFieldSplitSchurPreType*,Mat*);
extern PetscErrorCode PCFieldSplitSetSchurFactType(PC,PCFieldSplitSchurFactType);
extern PetscErrorCode PCFieldSplitGetSchurBlocks(PC,Mat*,Mat*,Mat*,Mat*);
extern PetscErrorCode PCFieldSplitSchurGetS(PC,Mat *S);
extern PetscErrorCode PCFieldSplitSchurRestoreS(PC,Mat *S);

extern PetscErrorCode PCGalerkinSetRestriction(PC,Mat);
extern PetscErrorCode PCGalerkinSetInterpolation(PC,Mat);

extern PetscErrorCode PCSetCoordinates(PC,PetscInt,PetscInt,PetscReal*);

extern PetscErrorCode PCPythonSetType(PC,const char[]);

extern PetscErrorCode PCSetDM(PC,DM);
extern PetscErrorCode PCGetDM(PC,DM*);

extern PetscErrorCode PCSetApplicationContext(PC,void*);
extern PetscErrorCode PCGetApplicationContext(PC,void*);

extern PetscErrorCode PCBiCGStabCUSPSetTolerance(PC,PetscReal);
extern PetscErrorCode PCBiCGStabCUSPSetIterations(PC,PetscInt);
extern PetscErrorCode PCBiCGStabCUSPSetUseVerboseMonitor(PC,PetscBool);

extern PetscErrorCode PCAINVCUSPSetDropTolerance(PC,PetscReal);
extern PetscErrorCode PCAINVCUSPUseScaling(PC,PetscBool);
extern PetscErrorCode PCAINVCUSPSetNonzeros(PC,PetscInt);
extern PetscErrorCode PCAINVCUSPSetLinParameter(PC,PetscInt);

extern PetscErrorCode PCPARMSSetGlobal(PC,PCPARMSGlobalType);
extern PetscErrorCode PCPARMSSetLocal(PC,PCPARMSLocalType);
extern PetscErrorCode PCPARMSSetSolveTolerances(PC,PetscReal,PetscInt);
extern PetscErrorCode PCPARMSSetSolveRestart(PC,PetscInt);
extern PetscErrorCode PCPARMSSetNonsymPerm(PC,PetscBool);
extern PetscErrorCode PCPARMSSetFill(PC,PetscInt,PetscInt,PetscInt);

extern PetscErrorCode PCGAMGSetType( PC,PCGAMGType);
extern PetscErrorCode PCGAMGGetType( PC,PCGAMGType*);
extern PetscErrorCode PCGAMGSetProcEqLim(PC,PetscInt);
extern PetscErrorCode PCGAMGSetRepartitioning(PC,PetscBool);
extern PetscErrorCode PCGAMGSetUseASMAggs(PC,PetscBool);
extern PetscErrorCode PCGAMGSetSolverType(PC,char[],PetscInt);
extern PetscErrorCode PCGAMGSetThreshold(PC,PetscReal);
extern PetscErrorCode PCGAMGSetCoarseEqLim(PC,PetscInt);
extern PetscErrorCode PCGAMGSetNlevels(PC,PetscInt);
extern PetscErrorCode PCGAMGSetNSmooths(PC,PetscInt);
extern PetscErrorCode PCGAMGSetSymGraph(PC,PetscBool);
extern PetscErrorCode PCGAMGSetSquareGraph(PC,PetscInt);
extern PetscErrorCode PCGAMGSetReuseInterpolation(PC,PetscBool);
extern PetscErrorCode PCGAMGFinalizePackage(void);
extern PetscErrorCode PCGAMGInitializePackage(void);
extern PetscErrorCode PCGAMGRegister(PCGAMGType,PetscErrorCode (*)(PC));

extern PetscErrorCode PCGAMGClassicalSetType(PC,PCGAMGClassicalType);
extern PetscErrorCode PCGAMGClassicalGetType(PC,PCGAMGClassicalType*);

extern PetscErrorCode PCBDDCSetChangeOfBasisMat(PC,Mat);
extern PetscErrorCode PCBDDCSetPrimalVerticesLocalIS(PC,IS);
extern PetscErrorCode PCBDDCSetCoarseningRatio(PC,PetscInt);
extern PetscErrorCode PCBDDCSetLevels(PC,PetscInt);
extern PetscErrorCode PCBDDCSetNullSpace(PC,MatNullSpace);
extern PetscErrorCode PCBDDCSetDirichletBoundaries(PC,IS);
extern PetscErrorCode PCBDDCSetDirichletBoundariesLocal(PC,IS);
extern PetscErrorCode PCBDDCGetDirichletBoundaries(PC,IS*);
extern PetscErrorCode PCBDDCGetDirichletBoundariesLocal(PC,IS*);
extern PetscErrorCode PCBDDCSetNeumannBoundaries(PC,IS);
extern PetscErrorCode PCBDDCSetNeumannBoundariesLocal(PC,IS);
extern PetscErrorCode PCBDDCGetNeumannBoundaries(PC,IS*);
extern PetscErrorCode PCBDDCGetNeumannBoundariesLocal(PC,IS*);
extern PetscErrorCode PCBDDCSetDofsSplitting(PC,PetscInt,IS[]);
extern PetscErrorCode PCBDDCSetDofsSplittingLocal(PC,PetscInt,IS[]);
extern PetscErrorCode PCBDDCSetLocalAdjacencyGraph(PC,PetscInt,const PetscInt[],const PetscInt[],PetscCopyMode);
extern PetscErrorCode PCBDDCCreateFETIDPOperators(PC,Mat*,PC*);
extern PetscErrorCode PCBDDCMatFETIDPGetRHS(Mat,Vec,Vec);
extern PetscErrorCode PCBDDCMatFETIDPGetSolution(Mat,Vec,Vec);

extern PetscErrorCode PCISSetUseStiffnessScaling(PC,PetscBool);
extern PetscErrorCode PCISSetSubdomainScalingFactor(PC,PetscScalar);
extern PetscErrorCode PCISSetSubdomainDiagonalScaling(PC,Vec);

extern PetscErrorCode PCMGSetType(PC,PCMGType);
extern PetscErrorCode PCMGGetType(PC,PCMGType*);
extern PetscErrorCode PCMGSetLevels(PC,PetscInt,MPI_Comm*);
extern PetscErrorCode PCMGGetLevels(PC,PetscInt*);

extern PetscErrorCode PCMGSetNumberSmoothUp(PC,PetscInt);
extern PetscErrorCode PCMGSetNumberSmoothDown(PC,PetscInt);
extern PetscErrorCode PCMGSetCycleType(PC,PCMGCycleType);
extern PetscErrorCode PCMGSetCycleTypeOnLevel(PC,PetscInt,PCMGCycleType);
extern PetscErrorCode PCMGSetCyclesOnLevel(PC,PetscInt,PetscInt);
extern PetscErrorCode PCMGMultiplicativeSetCycles(PC,PetscInt);
extern PetscErrorCode PCMGSetGalerkin(PC,PetscBool);
extern PetscErrorCode PCMGGetGalerkin(PC,PetscBool*);

extern PetscErrorCode PCMGSetRhs(PC,PetscInt,Vec);
extern PetscErrorCode PCMGSetX(PC,PetscInt,Vec);
extern PetscErrorCode PCMGSetR(PC,PetscInt,Vec);

extern PetscErrorCode PCMGSetRestriction(PC,PetscInt,Mat);
extern PetscErrorCode PCMGGetRestriction(PC,PetscInt,Mat*);
extern PetscErrorCode PCMGSetInterpolation(PC,PetscInt,Mat);
extern PetscErrorCode PCMGGetInterpolation(PC,PetscInt,Mat*);
extern PetscErrorCode PCMGSetRScale(PC,PetscInt,Vec);
extern PetscErrorCode PCMGGetRScale(PC,PetscInt,Vec*);
extern PetscErrorCode PCMGSetResidual(PC,PetscInt,PetscErrorCode (*)(Mat,Vec,Vec,Vec),Mat);
extern PetscErrorCode PCMGResidualDefault(Mat,Vec,Vec,Vec);
# 7 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscksp.h" 2

extern PetscErrorCode KSPInitializePackage(void);
# 23 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscksp.h"
typedef struct _p_KSP* KSP;
# 32 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscksp.h"
typedef const char* KSPType;
# 69 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscksp.h"
extern PetscClassId KSP_CLASSID;
extern PetscClassId DMKSP_CLASSID;

extern PetscErrorCode KSPCreate(MPI_Comm,KSP *);
extern PetscErrorCode KSPSetType(KSP,KSPType);
extern PetscErrorCode KSPGetType(KSP,KSPType *);
extern PetscErrorCode KSPSetUp(KSP);
extern PetscErrorCode KSPSetUpOnBlocks(KSP);
extern PetscErrorCode KSPSolve(KSP,Vec,Vec);
extern PetscErrorCode KSPSolveTranspose(KSP,Vec,Vec);
extern PetscErrorCode KSPReset(KSP);
extern PetscErrorCode KSPDestroy(KSP*);
extern PetscErrorCode KSPSetReusePreconditioner(KSP,PetscBool);
extern PetscErrorCode KSPSetSkipPCSetFromOptions(KSP,PetscBool);

extern PetscFunctionList KSPList;
extern PetscErrorCode KSPRegister(const char[],PetscErrorCode (*)(KSP));

extern PetscErrorCode KSPSetPCSide(KSP,PCSide);
extern PetscErrorCode KSPGetPCSide(KSP,PCSide*);
extern PetscErrorCode KSPSetTolerances(KSP,PetscReal,PetscReal,PetscReal,PetscInt);
extern PetscErrorCode KSPGetTolerances(KSP,PetscReal*,PetscReal*,PetscReal*,PetscInt*);
extern PetscErrorCode KSPSetInitialGuessNonzero(KSP,PetscBool );
extern PetscErrorCode KSPGetInitialGuessNonzero(KSP,PetscBool *);
extern PetscErrorCode KSPSetInitialGuessKnoll(KSP,PetscBool );
extern PetscErrorCode KSPGetInitialGuessKnoll(KSP,PetscBool *);
extern PetscErrorCode KSPSetErrorIfNotConverged(KSP,PetscBool );
extern PetscErrorCode KSPGetErrorIfNotConverged(KSP,PetscBool *);
extern PetscErrorCode KSPSetComputeEigenvalues(KSP,PetscBool );
extern PetscErrorCode KSPGetComputeEigenvalues(KSP,PetscBool *);
extern PetscErrorCode KSPSetComputeSingularValues(KSP,PetscBool );
extern PetscErrorCode KSPGetComputeSingularValues(KSP,PetscBool *);
extern PetscErrorCode KSPGetRhs(KSP,Vec *);
extern PetscErrorCode KSPGetSolution(KSP,Vec *);
extern PetscErrorCode KSPGetResidualNorm(KSP,PetscReal*);
extern PetscErrorCode KSPGetIterationNumber(KSP,PetscInt*);
extern PetscErrorCode KSPGetTotalIterations(KSP,PetscInt*);
extern PetscErrorCode KSPCreateVecs(KSP,PetscInt,Vec**,PetscInt,Vec**);
__attribute((deprecated)) static inline PetscErrorCode KSPGetVecs(KSP ksp,PetscInt n,Vec **x,PetscInt m,Vec **y) {return KSPCreateVecs(ksp,n,x,m,y);}

extern PetscErrorCode KSPSetPreSolve(KSP,PetscErrorCode (*)(KSP,Vec,Vec,void*),void*);
extern PetscErrorCode KSPSetPostSolve(KSP,PetscErrorCode (*)(KSP,Vec,Vec,void*),void*);

extern PetscErrorCode KSPSetPC(KSP,PC);
extern PetscErrorCode KSPGetPC(KSP,PC*);

extern PetscErrorCode KSPMonitor(KSP,PetscInt,PetscReal);
extern PetscErrorCode KSPMonitorSet(KSP,PetscErrorCode (*)(KSP,PetscInt,PetscReal,void*),void *,PetscErrorCode (*)(void**));
extern PetscErrorCode KSPMonitorCancel(KSP);
extern PetscErrorCode KSPGetMonitorContext(KSP,void **);
extern PetscErrorCode KSPGetResidualHistory(KSP,PetscReal*[],PetscInt *);
extern PetscErrorCode KSPSetResidualHistory(KSP,PetscReal[],PetscInt,PetscBool );

extern PetscErrorCode KSPBuildSolutionDefault(KSP,Vec,Vec*);
extern PetscErrorCode KSPBuildResidualDefault(KSP,Vec,Vec,Vec *);
extern PetscErrorCode KSPDestroyDefault(KSP);
extern PetscErrorCode KSPSetWorkVecs(KSP,PetscInt);

extern PetscErrorCode PCKSPGetKSP(PC,KSP*);
extern PetscErrorCode PCBJacobiGetSubKSP(PC,PetscInt*,PetscInt*,KSP*[]);
extern PetscErrorCode PCASMGetSubKSP(PC,PetscInt*,PetscInt*,KSP*[]);
extern PetscErrorCode PCGASMGetSubKSP(PC,PetscInt*,PetscInt*,KSP*[]);
extern PetscErrorCode PCFieldSplitGetSubKSP(PC,PetscInt*,KSP*[]);
extern PetscErrorCode PCMGGetSmoother(PC,PetscInt,KSP*);
extern PetscErrorCode PCMGGetSmootherDown(PC,PetscInt,KSP*);
extern PetscErrorCode PCMGGetSmootherUp(PC,PetscInt,KSP*);
extern PetscErrorCode PCMGGetCoarseSolve(PC,KSP*);
extern PetscErrorCode PCGalerkinGetKSP(PC,KSP *);

extern PetscErrorCode KSPBuildSolution(KSP,Vec,Vec *);
extern PetscErrorCode KSPBuildResidual(KSP,Vec,Vec,Vec *);

extern PetscErrorCode KSPRichardsonSetScale(KSP,PetscReal);
extern PetscErrorCode KSPRichardsonSetSelfScale(KSP,PetscBool );
extern PetscErrorCode KSPChebyshevSetEigenvalues(KSP,PetscReal,PetscReal);
extern PetscErrorCode KSPChebyshevEstEigSet(KSP,PetscReal,PetscReal,PetscReal,PetscReal);
extern PetscErrorCode KSPChebyshevEstEigSetRandom(KSP,PetscRandom);
extern PetscErrorCode KSPChebyshevEstEigGetKSP(KSP,KSP*);
extern PetscErrorCode KSPComputeExtremeSingularValues(KSP,PetscReal*,PetscReal*);
extern PetscErrorCode KSPComputeEigenvalues(KSP,PetscInt,PetscReal[],PetscReal[],PetscInt *);
extern PetscErrorCode KSPComputeEigenvaluesExplicitly(KSP,PetscInt,PetscReal[],PetscReal[]);
# 163 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscksp.h"
typedef enum {KSP_FCG_TRUNC_TYPE_STANDARD,KSP_FCG_TRUNC_TYPE_NOTAY} KSPFCGTruncationType;
extern const char *const KSPFCGTruncationTypes[];

extern PetscErrorCode KSPFCGSetMmax(KSP,PetscInt);
extern PetscErrorCode KSPFCGGetMmax(KSP,PetscInt*);
extern PetscErrorCode KSPFCGSetNprealloc(KSP,PetscInt);
extern PetscErrorCode KSPFCGGetNprealloc(KSP,PetscInt*);
extern PetscErrorCode KSPFCGSetTruncationType(KSP,KSPFCGTruncationType);
extern PetscErrorCode KSPFCGGetTruncationType(KSP,KSPFCGTruncationType*);

extern PetscErrorCode KSPGMRESSetRestart(KSP, PetscInt);
extern PetscErrorCode KSPGMRESGetRestart(KSP, PetscInt*);
extern PetscErrorCode KSPGMRESSetHapTol(KSP,PetscReal);

extern PetscErrorCode KSPGMRESSetPreAllocateVectors(KSP);
extern PetscErrorCode KSPGMRESSetOrthogonalization(KSP,PetscErrorCode (*)(KSP,PetscInt));
extern PetscErrorCode KSPGMRESGetOrthogonalization(KSP,PetscErrorCode (**)(KSP,PetscInt));
extern PetscErrorCode KSPGMRESModifiedGramSchmidtOrthogonalization(KSP,PetscInt);
extern PetscErrorCode KSPGMRESClassicalGramSchmidtOrthogonalization(KSP,PetscInt);

extern PetscErrorCode KSPLGMRESSetAugDim(KSP,PetscInt);
extern PetscErrorCode KSPLGMRESSetConstant(KSP);

extern PetscErrorCode KSPGCRSetRestart(KSP,PetscInt);
extern PetscErrorCode KSPGCRGetRestart(KSP,PetscInt*);
extern PetscErrorCode KSPGCRSetModifyPC(KSP,PetscErrorCode (*)(KSP,PetscInt,PetscReal,void*),void*,PetscErrorCode(*)(void*));
# 199 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscksp.h"
typedef enum {KSP_GMRES_CGS_REFINE_NEVER, KSP_GMRES_CGS_REFINE_IFNEEDED, KSP_GMRES_CGS_REFINE_ALWAYS} KSPGMRESCGSRefinementType;
extern const char *const KSPGMRESCGSRefinementTypes[];
# 243 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscksp.h"
extern PetscErrorCode KSPGMRESSetCGSRefinementType(KSP,KSPGMRESCGSRefinementType);
extern PetscErrorCode KSPGMRESGetCGSRefinementType(KSP,KSPGMRESCGSRefinementType*);

extern PetscErrorCode KSPFGMRESModifyPCNoChange(KSP,PetscInt,PetscInt,PetscReal,void*);
extern PetscErrorCode KSPFGMRESModifyPCKSP(KSP,PetscInt,PetscInt,PetscReal,void*);
extern PetscErrorCode KSPFGMRESSetModifyPC(KSP,PetscErrorCode (*)(KSP,PetscInt,PetscInt,PetscReal,void*),void*,PetscErrorCode(*)(void*));

extern PetscErrorCode KSPQCGSetTrustRegionRadius(KSP,PetscReal);
extern PetscErrorCode KSPQCGGetQuadratic(KSP,PetscReal*);
extern PetscErrorCode KSPQCGGetTrialStepNorm(KSP,PetscReal*);

extern PetscErrorCode KSPBCGSLSetXRes(KSP,PetscReal);
extern PetscErrorCode KSPBCGSLSetPol(KSP,PetscBool );
extern PetscErrorCode KSPBCGSLSetEll(KSP,PetscInt);
extern PetscErrorCode KSPBCGSLSetUsePseudoinverse(KSP,PetscBool);

extern PetscErrorCode KSPSetFromOptions(KSP);
extern PetscErrorCode KSPAddOptionsChecker(PetscErrorCode (*)(KSP));

extern PetscErrorCode KSPMonitorSingularValue(KSP,PetscInt,PetscReal,void *);
extern PetscErrorCode KSPMonitorDefault(KSP,PetscInt,PetscReal,void *);
extern PetscErrorCode KSPLSQRMonitorDefault(KSP,PetscInt,PetscReal,void *);
extern PetscErrorCode KSPMonitorRange(KSP,PetscInt,PetscReal,void *);
extern PetscErrorCode KSPMonitorDynamicTolerance(KSP ksp,PetscInt its,PetscReal fnorm,void *dummy);
extern PetscErrorCode KSPMonitorDynamicToleranceDestroy(void **dummy);
extern PetscErrorCode KSPMonitorTrueResidualNorm(KSP,PetscInt,PetscReal,void *);
extern PetscErrorCode KSPMonitorTrueResidualMaxNorm(KSP,PetscInt,PetscReal,void *);
extern PetscErrorCode KSPMonitorDefaultShort(KSP,PetscInt,PetscReal,void *);
extern PetscErrorCode KSPMonitorSolution(KSP,PetscInt,PetscReal,void *);
extern PetscErrorCode KSPMonitorSAWs(KSP,PetscInt,PetscReal,void*);
extern PetscErrorCode KSPMonitorSAWsCreate(KSP,void**);
extern PetscErrorCode KSPMonitorSAWsDestroy(void**);
extern PetscErrorCode KSPGMRESMonitorKrylov(KSP,PetscInt,PetscReal,void *);

extern PetscErrorCode KSPUnwindPreconditioner(KSP,Vec,Vec);
extern PetscErrorCode KSPInitialResidual(KSP,Vec,Vec,Vec,Vec,Vec);

extern PetscErrorCode KSPSetOperators(KSP,Mat,Mat);
extern PetscErrorCode KSPGetOperators(KSP,Mat*,Mat*);
extern PetscErrorCode KSPGetOperatorsSet(KSP,PetscBool *,PetscBool *);
extern PetscErrorCode KSPSetOptionsPrefix(KSP,const char[]);
extern PetscErrorCode KSPAppendOptionsPrefix(KSP,const char[]);
extern PetscErrorCode KSPGetOptionsPrefix(KSP,const char*[]);
extern PetscErrorCode KSPSetTabLevel(KSP,PetscInt);
extern PetscErrorCode KSPGetTabLevel(KSP,PetscInt*);

extern PetscErrorCode KSPSetDiagonalScale(KSP,PetscBool );
extern PetscErrorCode KSPGetDiagonalScale(KSP,PetscBool *);
extern PetscErrorCode KSPSetDiagonalScaleFix(KSP,PetscBool );
extern PetscErrorCode KSPGetDiagonalScaleFix(KSP,PetscBool *);

extern PetscErrorCode KSPView(KSP,PetscViewer);
extern PetscErrorCode KSPLoad(KSP,PetscViewer);
static inline PetscErrorCode KSPViewFromOptions(KSP A,PetscObject obj,const char name[]) {return PetscObjectViewFromOptions((PetscObject)A,obj,name);}
extern PetscErrorCode KSPReasonView(KSP,PetscViewer);
extern PetscErrorCode KSPReasonViewFromOptions(KSP);



extern PetscErrorCode KSPLSQRSetStandardErrorVec(KSP,Vec);
extern PetscErrorCode KSPLSQRGetStandardErrorVec(KSP,Vec*);

extern PetscErrorCode PCRedundantGetKSP(PC,KSP*);
extern PetscErrorCode PCRedistributeGetKSP(PC,KSP*);
# 322 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscksp.h"
typedef enum {KSP_NORM_DEFAULT = -1,KSP_NORM_NONE = 0,KSP_NORM_PRECONDITIONED = 1,KSP_NORM_UNPRECONDITIONED = 2,KSP_NORM_NATURAL = 3} KSPNormType;

extern const char *const*const KSPNormTypes;
# 365 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscksp.h"
extern PetscErrorCode KSPSetNormType(KSP,KSPNormType);
extern PetscErrorCode KSPGetNormType(KSP,KSPNormType*);
extern PetscErrorCode KSPSetSupportedNorm(KSP ksp,KSPNormType,PCSide,PetscInt);
extern PetscErrorCode KSPSetCheckNormIteration(KSP,PetscInt);
extern PetscErrorCode KSPSetLagNorm(KSP,PetscBool);
# 385 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscksp.h"
typedef enum {
              KSP_CONVERGED_RTOL_NORMAL = 1,
              KSP_CONVERGED_ATOL_NORMAL = 9,
              KSP_CONVERGED_RTOL = 2,
              KSP_CONVERGED_ATOL = 3,
              KSP_CONVERGED_ITS = 4,
              KSP_CONVERGED_CG_NEG_CURVE = 5,
              KSP_CONVERGED_CG_CONSTRAINED = 6,
              KSP_CONVERGED_STEP_LENGTH = 7,
              KSP_CONVERGED_HAPPY_BREAKDOWN = 8,

              KSP_DIVERGED_NULL = -2,
              KSP_DIVERGED_ITS = -3,
              KSP_DIVERGED_DTOL = -4,
              KSP_DIVERGED_BREAKDOWN = -5,
              KSP_DIVERGED_BREAKDOWN_BICG = -6,
              KSP_DIVERGED_NONSYMMETRIC = -7,
              KSP_DIVERGED_INDEFINITE_PC = -8,
              KSP_DIVERGED_NANORINF = -9,
              KSP_DIVERGED_INDEFINITE_MAT = -10,
              KSP_DIVERGED_PCSETUP_FAILED = -11,

              KSP_CONVERGED_ITERATING = 0} KSPConvergedReason;
extern const char *const*KSPConvergedReasons;
# 529 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscksp.h"
extern PetscErrorCode KSPSetConvergenceTest(KSP,PetscErrorCode (*)(KSP,PetscInt,PetscReal,KSPConvergedReason*,void*),void *,PetscErrorCode (*)(void*));
extern PetscErrorCode KSPGetConvergenceContext(KSP,void **);
extern PetscErrorCode KSPConvergedDefault(KSP,PetscInt,PetscReal,KSPConvergedReason*,void *);
extern PetscErrorCode KSPConvergedLSQR(KSP,PetscInt,PetscReal,KSPConvergedReason*,void *);
extern PetscErrorCode KSPConvergedDefaultDestroy(void *);
extern PetscErrorCode KSPConvergedDefaultCreate(void **);
extern PetscErrorCode KSPConvergedDefaultSetUIRNorm(KSP);
extern PetscErrorCode KSPConvergedDefaultSetUMIRNorm(KSP);
extern PetscErrorCode KSPConvergedSkip(KSP,PetscInt,PetscReal,KSPConvergedReason*,void *);
extern PetscErrorCode KSPGetConvergedReason(KSP,KSPConvergedReason *);

__attribute((deprecated)) static inline void KSPDefaultConverged(void) { }

__attribute((deprecated)) static inline void KSPDefaultConvergedDestroy(void) { }

__attribute((deprecated)) static inline void KSPDefaultConvergedCreate(void) { }

__attribute((deprecated)) static inline void KSPDefaultConvergedSetUIRNorm(void) { }

__attribute((deprecated)) static inline void KSPDefaultConvergedSetUMIRNorm(void) { }

__attribute((deprecated)) static inline void KSPSkipConverged(void) { }


extern PetscErrorCode KSPComputeExplicitOperator(KSP,Mat *);
# 562 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscksp.h"
typedef enum {KSP_CG_SYMMETRIC=0,KSP_CG_HERMITIAN=1} KSPCGType;
extern const char *const KSPCGTypes[];

extern PetscErrorCode KSPCGSetType(KSP,KSPCGType);
extern PetscErrorCode KSPCGUseSingleReduction(KSP,PetscBool );

extern PetscErrorCode KSPNASHSetRadius(KSP,PetscReal);
extern PetscErrorCode KSPNASHGetNormD(KSP,PetscReal *);
extern PetscErrorCode KSPNASHGetObjFcn(KSP,PetscReal *);

extern PetscErrorCode KSPSTCGSetRadius(KSP,PetscReal);
extern PetscErrorCode KSPSTCGGetNormD(KSP,PetscReal *);
extern PetscErrorCode KSPSTCGGetObjFcn(KSP,PetscReal *);

extern PetscErrorCode KSPGLTRSetRadius(KSP,PetscReal);
extern PetscErrorCode KSPGLTRGetNormD(KSP,PetscReal *);
extern PetscErrorCode KSPGLTRGetObjFcn(KSP,PetscReal *);
extern PetscErrorCode KSPGLTRGetMinEig(KSP,PetscReal *);
extern PetscErrorCode KSPGLTRGetLambda(KSP,PetscReal *);

extern PetscErrorCode KSPPythonSetType(KSP,const char[]);

extern PetscErrorCode PCPreSolve(PC,KSP);
extern PetscErrorCode PCPostSolve(PC,KSP);


extern PetscErrorCode KSPMonitorLGResidualNormCreate(const char[],const char[],int,int,int,int,PetscObject**);
extern PetscErrorCode KSPMonitorLGResidualNorm(KSP,PetscInt,PetscReal,PetscObject*);
extern PetscErrorCode KSPMonitorLGResidualNormDestroy(PetscObject**);
extern PetscErrorCode KSPMonitorLGTrueResidualNormCreate(const char[],const char[],int,int,int,int,PetscObject**);
extern PetscErrorCode KSPMonitorLGTrueResidualNorm(KSP,PetscInt,PetscReal,PetscObject*);
extern PetscErrorCode KSPMonitorLGTrueResidualNormDestroy(PetscObject**);
extern PetscErrorCode KSPMonitorLGRange(KSP,PetscInt,PetscReal,void*);

extern PetscErrorCode PCShellSetPreSolve(PC,PetscErrorCode (*)(PC,KSP,Vec,Vec));
extern PetscErrorCode PCShellSetPostSolve(PC,PetscErrorCode (*)(PC,KSP,Vec,Vec));


typedef struct _p_KSPFischerGuess {PetscInt method,curl,maxl,refcnt;PetscBool monitor;Mat mat; KSP ksp;}* KSPFischerGuess;

extern PetscErrorCode KSPFischerGuessCreate(KSP,PetscInt,PetscInt,KSPFischerGuess*);
extern PetscErrorCode KSPFischerGuessDestroy(KSPFischerGuess*);
extern PetscErrorCode KSPFischerGuessReset(KSPFischerGuess);
extern PetscErrorCode KSPFischerGuessUpdate(KSPFischerGuess,Vec);
extern PetscErrorCode KSPFischerGuessFormGuess(KSPFischerGuess,Vec,Vec);
extern PetscErrorCode KSPFischerGuessSetFromOptions(KSPFischerGuess);

extern PetscErrorCode KSPSetUseFischerGuess(KSP,PetscInt,PetscInt);
extern PetscErrorCode KSPSetFischerGuess(KSP,KSPFischerGuess);
extern PetscErrorCode KSPGetFischerGuess(KSP,KSPFischerGuess*);
# 620 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscksp.h"
typedef enum {MAT_SCHUR_COMPLEMENT_AINV_DIAG, MAT_SCHUR_COMPLEMENT_AINV_LUMP} MatSchurComplementAinvType;
extern const char *const MatSchurComplementAinvTypes[];

extern PetscErrorCode MatCreateSchurComplement(Mat,Mat,Mat,Mat,Mat,Mat*);
extern PetscErrorCode MatSchurComplementGetKSP(Mat,KSP*);
extern PetscErrorCode MatSchurComplementSetKSP(Mat,KSP);
extern PetscErrorCode MatSchurComplementSetSubMatrices(Mat,Mat,Mat,Mat,Mat,Mat);
extern PetscErrorCode MatSchurComplementUpdateSubMatrices(Mat,Mat,Mat,Mat,Mat,Mat);
extern PetscErrorCode MatSchurComplementGetSubMatrices(Mat,Mat*,Mat*,Mat*,Mat*,Mat*);
extern PetscErrorCode MatSchurComplementSetAinvType(Mat,MatSchurComplementAinvType);
extern PetscErrorCode MatSchurComplementGetAinvType(Mat,MatSchurComplementAinvType*);
extern PetscErrorCode MatSchurComplementGetPmat(Mat,MatReuse,Mat*);
extern PetscErrorCode MatSchurComplementComputeExplicitOperator(Mat,Mat*);
extern PetscErrorCode MatGetSchurComplement(Mat,IS,IS,IS,IS,MatReuse,Mat *,MatSchurComplementAinvType,MatReuse,Mat *);
extern PetscErrorCode MatCreateSchurComplementPmat(Mat,Mat,Mat,Mat,MatSchurComplementAinvType,MatReuse,Mat*);

extern PetscErrorCode KSPSetDM(KSP,DM);
extern PetscErrorCode KSPSetDMActive(KSP,PetscBool );
extern PetscErrorCode KSPGetDM(KSP,DM*);
extern PetscErrorCode KSPSetApplicationContext(KSP,void*);
extern PetscErrorCode KSPGetApplicationContext(KSP,void*);
extern PetscErrorCode KSPSetComputeRHS(KSP,PetscErrorCode (*func)(KSP,Vec,void*),void *);
extern PetscErrorCode KSPSetComputeOperators(KSP,PetscErrorCode(*)(KSP,Mat,Mat,void*),void*);
extern PetscErrorCode KSPSetComputeInitialGuess(KSP,PetscErrorCode(*)(KSP,Vec,void*),void*);
extern PetscErrorCode DMKSPSetComputeOperators(DM,PetscErrorCode(*)(KSP,Mat,Mat,void*),void*);
extern PetscErrorCode DMKSPGetComputeOperators(DM,PetscErrorCode(**)(KSP,Mat,Mat,void*),void*);
extern PetscErrorCode DMKSPSetComputeRHS(DM,PetscErrorCode(*)(KSP,Vec,void*),void*);
extern PetscErrorCode DMKSPGetComputeRHS(DM,PetscErrorCode(**)(KSP,Vec,void*),void*);
extern PetscErrorCode DMKSPSetComputeInitialGuess(DM,PetscErrorCode(*)(KSP,Vec,void*),void*);
extern PetscErrorCode DMKSPGetComputeInitialGuess(DM,PetscErrorCode(**)(KSP,Vec,void*),void*);

extern PetscErrorCode DMGlobalToLocalSolve(DM,Vec,Vec);
extern PetscErrorCode DMPlexProjectField(DM, Vec, void (**)(const PetscScalar[], const PetscScalar[], const PetscScalar[], const PetscScalar[], const PetscScalar[], const PetscScalar[], const PetscReal [], PetscScalar []), InsertMode, Vec);
# 4 "variables.h" 2
# 1 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsnes.h" 1







# 1 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscfvtypes.h" 1
# 13 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscfvtypes.h"
typedef struct _p_PetscLimiter *PetscLimiter;
# 24 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscfvtypes.h"
typedef struct _p_PetscFV *PetscFV;
# 40 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscfvtypes.h"
typedef struct {
  PetscReal normal[3];
  PetscReal centroid[3];
  PetscScalar grad[2][3];
} PetscFVFaceGeom;
# 59 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscfvtypes.h"
typedef struct {
  PetscReal centroid[3];
  PetscReal volume;
} PetscFVCellGeom;
# 9 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsnes.h" 2
# 20 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsnes.h"
typedef struct _p_SNES* SNES;
# 29 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsnes.h"
typedef const char* SNESType;
# 51 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsnes.h"
extern PetscClassId SNES_CLASSID;
extern PetscClassId DMSNES_CLASSID;

extern PetscErrorCode SNESInitializePackage(void);

extern PetscErrorCode SNESCreate(MPI_Comm,SNES*);
extern PetscErrorCode SNESReset(SNES);
extern PetscErrorCode SNESDestroy(SNES*);
extern PetscErrorCode SNESSetType(SNES,SNESType);
extern PetscErrorCode SNESMonitor(SNES,PetscInt,PetscReal);
extern PetscErrorCode SNESMonitorSet(SNES,PetscErrorCode(*)(SNES,PetscInt,PetscReal,void*),void *,PetscErrorCode (*)(void**));
extern PetscErrorCode SNESMonitorCancel(SNES);
extern PetscErrorCode SNESMonitorSAWs(SNES,PetscInt,PetscReal,void*);
extern PetscErrorCode SNESMonitorSAWsCreate(SNES,void**);
extern PetscErrorCode SNESMonitorSAWsDestroy(void**);
extern PetscErrorCode SNESSetConvergenceHistory(SNES,PetscReal[],PetscInt[],PetscInt,PetscBool );
extern PetscErrorCode SNESGetConvergenceHistory(SNES,PetscReal*[],PetscInt *[],PetscInt *);
extern PetscErrorCode SNESSetUp(SNES);
extern PetscErrorCode SNESSolve(SNES,Vec,Vec);
extern PetscErrorCode SNESSetErrorIfNotConverged(SNES,PetscBool );
extern PetscErrorCode SNESGetErrorIfNotConverged(SNES,PetscBool *);

extern PetscErrorCode SNESSetWorkVecs(SNES,PetscInt);

extern PetscErrorCode SNESAddOptionsChecker(PetscErrorCode (*)(SNES));

extern PetscErrorCode SNESSetUpdate(SNES, PetscErrorCode (*)(SNES, PetscInt));


extern PetscErrorCode SNESRegister(const char[],PetscErrorCode (*)(SNES));

extern PetscErrorCode SNESGetKSP(SNES,KSP*);
extern PetscErrorCode SNESSetKSP(SNES,KSP);
extern PetscErrorCode SNESSetSolution(SNES,Vec);
extern PetscErrorCode SNESGetSolution(SNES,Vec*);
extern PetscErrorCode SNESGetSolutionUpdate(SNES,Vec*);
extern PetscErrorCode SNESGetRhs(SNES,Vec*);
extern PetscErrorCode SNESView(SNES,PetscViewer);
extern PetscErrorCode SNESLoad(SNES,PetscViewer);
static inline PetscErrorCode SNESViewFromOptions(SNES A,PetscObject obj,const char name[]) {return PetscObjectViewFromOptions((PetscObject)A,obj,name);}
extern PetscErrorCode SNESReasonView(SNES,PetscViewer);
extern PetscErrorCode SNESReasonViewFromOptions(SNES);



extern PetscErrorCode SNESSetOptionsPrefix(SNES,const char[]);
extern PetscErrorCode SNESAppendOptionsPrefix(SNES,const char[]);
extern PetscErrorCode SNESGetOptionsPrefix(SNES,const char*[]);
extern PetscErrorCode SNESSetFromOptions(SNES);

extern PetscErrorCode MatCreateSNESMF(SNES,Mat*);
extern PetscErrorCode MatMFFDComputeJacobian(SNES,Vec,Mat,Mat,void*);

extern PetscErrorCode MatDAADSetSNES(Mat,SNES);

extern PetscErrorCode SNESGetType(SNES,SNESType*);
extern PetscErrorCode SNESMonitorDefault(SNES,PetscInt,PetscReal,void *);
extern PetscErrorCode SNESMonitorRange(SNES,PetscInt,PetscReal,void *);
extern PetscErrorCode SNESMonitorRatio(SNES,PetscInt,PetscReal,void *);
extern PetscErrorCode SNESMonitorSetRatio(SNES,PetscViewer);
extern PetscErrorCode SNESMonitorSolution(SNES,PetscInt,PetscReal,void *);
extern PetscErrorCode SNESMonitorResidual(SNES,PetscInt,PetscReal,void *);
extern PetscErrorCode SNESMonitorSolutionUpdate(SNES,PetscInt,PetscReal,void *);
extern PetscErrorCode SNESMonitorDefaultShort(SNES,PetscInt,PetscReal,void *);
extern PetscErrorCode SNESMonitorDefaultField(SNES,PetscInt,PetscReal,void *);
extern PetscErrorCode SNESMonitorJacUpdateSpectrum(SNES,PetscInt,PetscReal,void *);
extern PetscErrorCode SNESMonitorFields(SNES,PetscInt,PetscReal,void *);
extern PetscErrorCode KSPMonitorSNES(KSP,PetscInt,PetscReal,void*);
extern PetscErrorCode KSPMonitorSNESLGResidualNormCreate(const char[],const char[],int,int,int,int,PetscObject**);
extern PetscErrorCode KSPMonitorSNESLGResidualNorm(KSP,PetscInt,PetscReal,PetscObject*);
extern PetscErrorCode KSPMonitorSNESLGResidualNormDestroy(PetscObject**);

extern PetscErrorCode SNESSetTolerances(SNES,PetscReal,PetscReal,PetscReal,PetscInt,PetscInt);
extern PetscErrorCode SNESGetTolerances(SNES,PetscReal*,PetscReal*,PetscReal*,PetscInt*,PetscInt*);
extern PetscErrorCode SNESSetTrustRegionTolerance(SNES,PetscReal);
extern PetscErrorCode SNESGetIterationNumber(SNES,PetscInt*);
extern PetscErrorCode SNESSetIterationNumber(SNES,PetscInt);

extern PetscErrorCode SNESGetNonlinearStepFailures(SNES,PetscInt*);
extern PetscErrorCode SNESSetMaxNonlinearStepFailures(SNES,PetscInt);
extern PetscErrorCode SNESGetMaxNonlinearStepFailures(SNES,PetscInt*);
extern PetscErrorCode SNESGetNumberFunctionEvals(SNES,PetscInt*);

extern PetscErrorCode SNESSetLagPreconditioner(SNES,PetscInt);
extern PetscErrorCode SNESGetLagPreconditioner(SNES,PetscInt*);
extern PetscErrorCode SNESSetLagJacobian(SNES,PetscInt);
extern PetscErrorCode SNESGetLagJacobian(SNES,PetscInt*);
extern PetscErrorCode SNESSetLagPreconditionerPersists(SNES,PetscBool);
extern PetscErrorCode SNESSetLagJacobianPersists(SNES,PetscBool);
extern PetscErrorCode SNESSetGridSequence(SNES,PetscInt);
extern PetscErrorCode SNESGetGridSequence(SNES,PetscInt*);

extern PetscErrorCode SNESGetLinearSolveIterations(SNES,PetscInt*);
extern PetscErrorCode SNESGetLinearSolveFailures(SNES,PetscInt*);
extern PetscErrorCode SNESSetMaxLinearSolveFailures(SNES,PetscInt);
extern PetscErrorCode SNESGetMaxLinearSolveFailures(SNES,PetscInt*);
extern PetscErrorCode SNESSetCountersReset(SNES,PetscBool);

extern PetscErrorCode SNESKSPSetUseEW(SNES,PetscBool );
extern PetscErrorCode SNESKSPGetUseEW(SNES,PetscBool *);
extern PetscErrorCode SNESKSPSetParametersEW(SNES,PetscInt,PetscReal,PetscReal,PetscReal,PetscReal,PetscReal,PetscReal);
extern PetscErrorCode SNESKSPGetParametersEW(SNES,PetscInt*,PetscReal*,PetscReal*,PetscReal*,PetscReal*,PetscReal*,PetscReal*);


extern PetscErrorCode SNESMonitorLGCreate(const char[],const char[],int,int,int,int,PetscObject**);
extern PetscErrorCode SNESMonitorLGResidualNorm(SNES,PetscInt,PetscReal,PetscObject*);
extern PetscErrorCode SNESMonitorLGDestroy(PetscObject**);
extern PetscErrorCode SNESMonitorLGRange(SNES,PetscInt,PetscReal,void*);

extern PetscErrorCode SNESSetApplicationContext(SNES,void *);
extern PetscErrorCode SNESGetApplicationContext(SNES,void *);
extern PetscErrorCode SNESSetComputeApplicationContext(SNES,PetscErrorCode (*)(SNES,void**),PetscErrorCode (*)(void**));

extern PetscErrorCode SNESPythonSetType(SNES,const char[]);

extern PetscErrorCode SNESSetFunctionDomainError(SNES);
extern PetscErrorCode SNESGetFunctionDomainError(SNES, PetscBool *);
# 218 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsnes.h"
typedef enum {
              SNES_CONVERGED_FNORM_ABS = 2,
              SNES_CONVERGED_FNORM_RELATIVE = 3,
              SNES_CONVERGED_SNORM_RELATIVE = 4,
              SNES_CONVERGED_ITS = 5,
              SNES_CONVERGED_TR_DELTA = 7,

              SNES_DIVERGED_FUNCTION_DOMAIN = -1,
              SNES_DIVERGED_FUNCTION_COUNT = -2,
              SNES_DIVERGED_LINEAR_SOLVE = -3,
              SNES_DIVERGED_FNORM_NAN = -4,
              SNES_DIVERGED_MAX_IT = -5,
              SNES_DIVERGED_LINE_SEARCH = -6,
              SNES_DIVERGED_INNER = -7,
              SNES_DIVERGED_LOCAL_MIN = -8,
              SNES_CONVERGED_ITERATING = 0} SNESConvergedReason;
extern const char *const*SNESConvergedReasons;
# 324 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsnes.h"
extern PetscErrorCode SNESSetConvergenceTest(SNES,PetscErrorCode (*)(SNES,PetscInt,PetscReal,PetscReal,PetscReal,SNESConvergedReason*,void*),void*,PetscErrorCode (*)(void*));
extern PetscErrorCode SNESConvergedDefault(SNES,PetscInt,PetscReal,PetscReal,PetscReal,SNESConvergedReason*,void*);
extern PetscErrorCode SNESConvergedSkip(SNES,PetscInt,PetscReal,PetscReal,PetscReal,SNESConvergedReason*,void*);
extern PetscErrorCode SNESGetConvergedReason(SNES,SNESConvergedReason*);

__attribute((deprecated)) static inline void SNESSkipConverged(void) { }



extern PetscErrorCode SNESSetFunction(SNES,Vec,PetscErrorCode (*)(SNES,Vec,Vec,void*),void*);
extern PetscErrorCode SNESGetFunction(SNES,Vec*,PetscErrorCode (**)(SNES,Vec,Vec,void*),void**);
extern PetscErrorCode SNESComputeFunction(SNES,Vec,Vec);

extern PetscErrorCode SNESSetJacobian(SNES,Mat,Mat,PetscErrorCode (*)(SNES,Vec,Mat,Mat,void*),void*);
extern PetscErrorCode SNESGetJacobian(SNES,Mat*,Mat*,PetscErrorCode (**)(SNES,Vec,Mat,Mat,void*),void**);
extern PetscErrorCode SNESObjectiveComputeFunctionDefaultFD(SNES,Vec,Vec,void*);
extern PetscErrorCode SNESComputeJacobianDefault(SNES,Vec,Mat,Mat,void*);
extern PetscErrorCode SNESComputeJacobianDefaultColor(SNES,Vec,Mat,Mat,void*);
extern PetscErrorCode SNESSetComputeInitialGuess(SNES,PetscErrorCode (*)(SNES,Vec,void*),void*);
extern PetscErrorCode SNESSetPicard(SNES,Vec,PetscErrorCode (*)(SNES,Vec,Vec,void*),Mat,Mat,PetscErrorCode (*)(SNES,Vec,Mat,Mat,void*),void*);
extern PetscErrorCode SNESGetPicard(SNES,Vec*,PetscErrorCode (**)(SNES,Vec,Vec,void*),Mat*,Mat*,PetscErrorCode (**)(SNES,Vec,Mat,Mat,void*),void**);
extern PetscErrorCode SNESSetInitialFunction(SNES,Vec);

extern PetscErrorCode SNESSetObjective(SNES,PetscErrorCode (*)(SNES,Vec,PetscReal *,void*),void*);
extern PetscErrorCode SNESGetObjective(SNES,PetscErrorCode (**)(SNES,Vec,PetscReal *,void*),void**);
extern PetscErrorCode SNESComputeObjective(SNES,Vec,PetscReal *);
# 366 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsnes.h"
typedef enum {SNES_NORM_DEFAULT = -1,
              SNES_NORM_NONE = 0,
              SNES_NORM_ALWAYS = 1,
              SNES_NORM_INITIAL_ONLY = 2,
              SNES_NORM_FINAL_ONLY = 3,
              SNES_NORM_INITIAL_FINAL_ONLY = 4} SNESNormSchedule;
extern const char *const*const SNESNormSchedules;
# 436 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsnes.h"
extern PetscErrorCode SNESSetNormSchedule(SNES,SNESNormSchedule);
extern PetscErrorCode SNESGetNormSchedule(SNES,SNESNormSchedule*);
# 449 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsnes.h"
typedef enum {SNES_FUNCTION_DEFAULT = -1,
              SNES_FUNCTION_UNPRECONDITIONED = 0,
              SNES_FUNCTION_PRECONDITIONED = 1} SNESFunctionType;
extern const char *const*const SNESFunctionTypes;

extern PetscErrorCode SNESSetFunctionType(SNES,SNESFunctionType);
extern PetscErrorCode SNESGetFunctionType(SNES,SNESFunctionType*);

extern PetscErrorCode SNESSetNGS(SNES,PetscErrorCode (*)(SNES,Vec,Vec,void*),void*);
extern PetscErrorCode SNESGetNGS(SNES,PetscErrorCode (**)(SNES,Vec,Vec,void*),void**);
extern PetscErrorCode SNESSetUseNGS(SNES,PetscBool);
extern PetscErrorCode SNESGetUseNGS(SNES,PetscBool *);
extern PetscErrorCode SNESComputeNGS(SNES,Vec,Vec);

extern PetscErrorCode SNESNGSSetSweeps(SNES,PetscInt);
extern PetscErrorCode SNESNGSGetSweeps(SNES,PetscInt *);
extern PetscErrorCode SNESNGSSetTolerances(SNES,PetscReal,PetscReal,PetscReal,PetscInt);
extern PetscErrorCode SNESNGSGetTolerances(SNES,PetscReal*,PetscReal*,PetscReal*,PetscInt*);

extern PetscErrorCode SNESUpdateCheckJacobian(SNES,PetscInt);

extern PetscErrorCode SNESShellGetContext(SNES,void**);
extern PetscErrorCode SNESShellSetContext(SNES,void*);
extern PetscErrorCode SNESShellSetSolve(SNES,PetscErrorCode (*)(SNES,Vec));
# 486 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsnes.h"
typedef struct _p_LineSearch* SNESLineSearch;
# 495 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsnes.h"
typedef const char* SNESLineSearchType;







extern PetscFunctionList SNESList;
extern PetscClassId SNESLINESEARCH_CLASSID;
extern PetscFunctionList SNESLineSearchList;
extern PetscLogEvent SNESLineSearch_Apply;





 typedef PetscErrorCode (*SNESLineSearchVIProjectFunc)(SNES,Vec);
 typedef PetscErrorCode (*SNESLineSearchVINormFunc)(SNES,Vec,Vec,PetscReal *);
 typedef PetscErrorCode (*SNESLineSearchApplyFunc)(SNESLineSearch);
 typedef PetscErrorCode (*SNESLineSearchUserFunc)(SNESLineSearch, void *);

extern PetscErrorCode SNESLineSearchCreate(MPI_Comm, SNESLineSearch*);
extern PetscErrorCode SNESLineSearchReset(SNESLineSearch);
extern PetscErrorCode SNESLineSearchView(SNESLineSearch,PetscViewer);
extern PetscErrorCode SNESLineSearchDestroy(SNESLineSearch *);
extern PetscErrorCode SNESLineSearchSetType(SNESLineSearch, SNESLineSearchType);
extern PetscErrorCode SNESLineSearchSetFromOptions(SNESLineSearch);
extern PetscErrorCode SNESLineSearchSetFunction(SNESLineSearch,PetscErrorCode (*)(SNES,Vec,Vec));
extern PetscErrorCode SNESLineSearchSetUp(SNESLineSearch);
extern PetscErrorCode SNESLineSearchApply(SNESLineSearch, Vec, Vec, PetscReal *, Vec);
extern PetscErrorCode SNESLineSearchPreCheck(SNESLineSearch,Vec,Vec,PetscBool *);
extern PetscErrorCode SNESLineSearchPostCheck(SNESLineSearch,Vec,Vec,Vec,PetscBool *,PetscBool *);
extern PetscErrorCode SNESLineSearchSetWorkVecs(SNESLineSearch, PetscInt);



extern PetscErrorCode SNESLineSearchSetPreCheck(SNESLineSearch, PetscErrorCode (*)(SNESLineSearch,Vec,Vec,PetscBool*,void*), void *ctx);
extern PetscErrorCode SNESLineSearchSetPostCheck(SNESLineSearch, PetscErrorCode (*)(SNESLineSearch,Vec,Vec,Vec,PetscBool *,PetscBool *,void*), void *ctx);

extern PetscErrorCode SNESLineSearchGetPreCheck(SNESLineSearch, PetscErrorCode (**)(SNESLineSearch,Vec,Vec,PetscBool*,void*), void **ctx);
extern PetscErrorCode SNESLineSearchGetPostCheck(SNESLineSearch, PetscErrorCode (**)(SNESLineSearch,Vec,Vec,Vec,PetscBool *,PetscBool *,void*), void **ctx);



extern PetscErrorCode SNESLineSearchSetVIFunctions(SNESLineSearch, SNESLineSearchVIProjectFunc, SNESLineSearchVINormFunc);
extern PetscErrorCode SNESLineSearchGetVIFunctions(SNESLineSearch, SNESLineSearchVIProjectFunc*, SNESLineSearchVINormFunc*);


extern PetscErrorCode SNESLineSearchSetSNES(SNESLineSearch,SNES);
extern PetscErrorCode SNESLineSearchGetSNES(SNESLineSearch,SNES*);


extern PetscErrorCode SNESLineSearchGetTolerances(SNESLineSearch,PetscReal*,PetscReal*,PetscReal*,PetscReal*,PetscReal*,PetscInt*);
extern PetscErrorCode SNESLineSearchSetTolerances(SNESLineSearch,PetscReal,PetscReal,PetscReal,PetscReal,PetscReal,PetscInt);

extern PetscErrorCode SNESLineSearchPreCheckPicard(SNESLineSearch,Vec,Vec,PetscBool*,void*);

extern PetscErrorCode SNESLineSearchGetLambda(SNESLineSearch,PetscReal*);
extern PetscErrorCode SNESLineSearchSetLambda(SNESLineSearch,PetscReal);

extern PetscErrorCode SNESLineSearchGetDamping(SNESLineSearch,PetscReal*);
extern PetscErrorCode SNESLineSearchSetDamping(SNESLineSearch,PetscReal);

extern PetscErrorCode SNESLineSearchGetOrder(SNESLineSearch,PetscInt *order);
extern PetscErrorCode SNESLineSearchSetOrder(SNESLineSearch,PetscInt order);
# 574 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsnes.h"
typedef enum {SNES_LINESEARCH_SUCCEEDED,
              SNES_LINESEARCH_FAILED_NANORINF,
              SNES_LINESEARCH_FAILED_DOMAIN,
              SNES_LINESEARCH_FAILED_REDUCT,
              SNES_LINESEARCH_FAILED_USER,
              SNES_LINESEARCH_FAILED_FUNCTION} SNESLineSearchReason;

extern PetscErrorCode SNESLineSearchGetReason(SNESLineSearch, SNESLineSearchReason*);
extern PetscErrorCode SNESLineSearchSetReason(SNESLineSearch, SNESLineSearchReason);

extern PetscErrorCode SNESLineSearchGetVecs(SNESLineSearch,Vec*,Vec*,Vec*,Vec*,Vec*);
extern PetscErrorCode SNESLineSearchSetVecs(SNESLineSearch,Vec,Vec,Vec,Vec,Vec);

extern PetscErrorCode SNESLineSearchGetNorms(SNESLineSearch, PetscReal *, PetscReal *, PetscReal *);
extern PetscErrorCode SNESLineSearchSetNorms(SNESLineSearch, PetscReal, PetscReal, PetscReal);
extern PetscErrorCode SNESLineSearchComputeNorms(SNESLineSearch);
extern PetscErrorCode SNESLineSearchSetComputeNorms(SNESLineSearch, PetscBool);

extern PetscErrorCode SNESLineSearchSetMonitor(SNESLineSearch, PetscBool);
extern PetscErrorCode SNESLineSearchGetMonitor(SNESLineSearch, PetscViewer*);

extern PetscErrorCode SNESLineSearchAppendOptionsPrefix(SNESLineSearch, const char prefix[]);
extern PetscErrorCode SNESLineSearchGetOptionsPrefix(SNESLineSearch, const char *prefix[]);



extern PetscErrorCode SNESLineSearchShellSetUserFunc(SNESLineSearch,SNESLineSearchUserFunc,void*);
extern PetscErrorCode SNESLineSearchShellGetUserFunc(SNESLineSearch,SNESLineSearchUserFunc*,void**);


extern PetscErrorCode SNESLineSearchBTSetAlpha(SNESLineSearch, PetscReal);
extern PetscErrorCode SNESLineSearchBTGetAlpha(SNESLineSearch, PetscReal*);


extern PetscErrorCode SNESLineSearchRegister(const char[],PetscErrorCode(*)(SNESLineSearch));


extern PetscErrorCode SNESVISetVariableBounds(SNES,Vec,Vec);
extern PetscErrorCode SNESVISetComputeVariableBounds(SNES, PetscErrorCode (*)(SNES,Vec,Vec));
extern PetscErrorCode SNESVIGetInactiveSet(SNES,IS*);
extern PetscErrorCode SNESVIGetActiveSetIS(SNES,Vec,Vec,IS*);
extern PetscErrorCode SNESVIComputeInactiveSetFnorm(SNES,Vec,Vec,PetscReal*);
extern PetscErrorCode SNESVISetRedundancyCheck(SNES,PetscErrorCode(*)(SNES,IS,IS*,void*),void*);

extern PetscErrorCode SNESTestLocalMin(SNES);


extern PetscErrorCode SNESComputeJacobian(SNES,Vec,Mat,Mat);

extern PetscErrorCode SNESSetDM(SNES,DM);
extern PetscErrorCode SNESGetDM(SNES,DM*);
extern PetscErrorCode SNESSetNPC(SNES,SNES);
extern PetscErrorCode SNESGetNPC(SNES,SNES*);
extern PetscErrorCode SNESHasNPC(SNES,PetscBool*);
extern PetscErrorCode SNESApplyNPC(SNES,Vec,Vec,Vec);
extern PetscErrorCode SNESGetNPCFunction(SNES,Vec,PetscReal*);
extern PetscErrorCode SNESComputeFunctionDefaultNPC(SNES,Vec,Vec);
extern PetscErrorCode SNESSetNPCSide(SNES,PCSide);
extern PetscErrorCode SNESGetNPCSide(SNES,PCSide*);
extern PetscErrorCode SNESSetLineSearch(SNES,SNESLineSearch);
extern PetscErrorCode SNESGetLineSearch(SNES,SNESLineSearch*);
extern PetscErrorCode SNESRestrictHookAdd(SNES,PetscErrorCode (*)(SNES,SNES,void*),void*);
extern PetscErrorCode SNESRestrictHooksRun(SNES,SNES);

__attribute((deprecated)) static inline PetscErrorCode SNESGetSNESLineSearch(SNES snes,SNESLineSearch *ls) {return SNESGetLineSearch(snes,ls);}
__attribute((deprecated)) static inline PetscErrorCode SNESSetSNESLineSearch(SNES snes,SNESLineSearch ls) {return SNESSetLineSearch(snes,ls);}

extern PetscErrorCode SNESSetUpMatrices(SNES);
extern PetscErrorCode DMSNESSetFunction(DM,PetscErrorCode(*)(SNES,Vec,Vec,void*),void*);
extern PetscErrorCode DMSNESGetFunction(DM,PetscErrorCode(**)(SNES,Vec,Vec,void*),void**);
extern PetscErrorCode DMSNESSetNGS(DM,PetscErrorCode(*)(SNES,Vec,Vec,void*),void*);
extern PetscErrorCode DMSNESGetNGS(DM,PetscErrorCode(**)(SNES,Vec,Vec,void*),void**);
extern PetscErrorCode DMSNESSetJacobian(DM,PetscErrorCode(*)(SNES,Vec,Mat,Mat,void*),void*);
extern PetscErrorCode DMSNESGetJacobian(DM,PetscErrorCode(**)(SNES,Vec,Mat,Mat,void*),void**);
extern PetscErrorCode DMSNESSetPicard(DM,PetscErrorCode(*)(SNES,Vec,Vec,void*),PetscErrorCode(*)(SNES,Vec,Mat,Mat,void*),void*);
extern PetscErrorCode DMSNESGetPicard(DM,PetscErrorCode(**)(SNES,Vec,Vec,void*),PetscErrorCode(**)(SNES,Vec,Mat,Mat,void*),void**);
extern PetscErrorCode DMSNESSetObjective(DM,PetscErrorCode (*)(SNES,Vec,PetscReal *,void*),void*);
extern PetscErrorCode DMSNESGetObjective(DM,PetscErrorCode (**)(SNES,Vec,PetscReal *,void*),void**);

 typedef PetscErrorCode (*DMDASNESFunction)(DMDALocalInfo*,void*,void*,void*);
 typedef PetscErrorCode (*DMDASNESJacobian)(DMDALocalInfo*,void*,Mat,Mat,void*);
 typedef PetscErrorCode (*DMDASNESObjective)(DMDALocalInfo*,void*,PetscReal*,void*);

extern PetscErrorCode DMDASNESSetFunctionLocal(DM,InsertMode,DMDASNESFunction,void*);
extern PetscErrorCode DMDASNESSetJacobianLocal(DM,DMDASNESJacobian,void*);
extern PetscErrorCode DMDASNESSetObjectiveLocal(DM,DMDASNESObjective,void*);
extern PetscErrorCode DMDASNESSetPicardLocal(DM,InsertMode,PetscErrorCode (*)(DMDALocalInfo*,void*,void*,void*),PetscErrorCode (*)(DMDALocalInfo*,void*,Mat,Mat,void*),void*);

extern PetscErrorCode DMPlexSNESGetGeometryFEM(DM,Vec*);
extern PetscErrorCode DMPlexSNESGetGeometryFVM(DM,Vec*,Vec*,PetscReal*);
extern PetscErrorCode DMPlexSNESGetGradientDM(DM,PetscFV,DM*);
extern PetscErrorCode DMPlexGetCellFields(DM, PetscInt, PetscInt, Vec, Vec, Vec, PetscScalar **, PetscScalar **, PetscScalar **);
extern PetscErrorCode DMPlexRestoreCellFields(DM, PetscInt, PetscInt, Vec, Vec, Vec, PetscScalar **, PetscScalar **, PetscScalar **);
extern PetscErrorCode DMPlexGetFaceFields(DM, PetscInt, PetscInt, Vec, Vec, Vec, Vec, Vec, PetscScalar **, PetscScalar **);
extern PetscErrorCode DMPlexRestoreFaceFields(DM, PetscInt, PetscInt, Vec, Vec, Vec, Vec, Vec, PetscScalar **, PetscScalar **);
extern PetscErrorCode DMPlexGetFaceGeometry(DM, PetscInt, PetscInt, Vec, Vec, PetscFVFaceGeom **, PetscReal **);
extern PetscErrorCode DMPlexRestoreFaceGeometry(DM, PetscInt, PetscInt, Vec, Vec, PetscFVFaceGeom **, PetscReal **);

extern PetscErrorCode DMSNESSetFunctionLocal(DM,PetscErrorCode (*)(DM,Vec,Vec,void*),void*);
extern PetscErrorCode DMSNESSetJacobianLocal(DM,PetscErrorCode (*)(DM,Vec,Mat,Mat,void*),void*);


extern PetscErrorCode SNESMultiblockSetFields(SNES, const char [], PetscInt, const PetscInt *);
extern PetscErrorCode SNESMultiblockSetIS(SNES, const char [], IS);
extern PetscErrorCode SNESMultiblockSetBlockSize(SNES, PetscInt);
extern PetscErrorCode SNESMultiblockSetType(SNES, PCCompositeType);
# 688 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsnes.h"
typedef const char* SNESMSType;
# 698 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsnes.h"
extern PetscErrorCode SNESMSRegister(SNESMSType,PetscInt,PetscInt,PetscReal,const PetscReal[],const PetscReal[],const PetscReal[]);
extern PetscErrorCode SNESMSSetType(SNES,SNESMSType);
extern PetscErrorCode SNESMSFinalizePackage(void);
extern PetscErrorCode SNESMSInitializePackage(void);
extern PetscErrorCode SNESMSRegisterDestroy(void);



typedef enum {
  SNES_NGMRES_RESTART_NONE = 0,
  SNES_NGMRES_RESTART_PERIODIC = 1,
  SNES_NGMRES_RESTART_DIFFERENCE = 2} SNESNGMRESRestartType;
extern const char *const SNESNGMRESRestartTypes[];

typedef enum {
  SNES_NGMRES_SELECT_NONE = 0,
  SNES_NGMRES_SELECT_DIFFERENCE = 1,
  SNES_NGMRES_SELECT_LINESEARCH = 2} SNESNGMRESSelectType;
extern const char *const SNESNGMRESSelectTypes[];

extern PetscErrorCode SNESNGMRESSetRestartType(SNES, SNESNGMRESRestartType);
extern PetscErrorCode SNESNGMRESSetSelectType(SNES, SNESNGMRESSelectType);



typedef enum {
  SNES_NCG_FR = 0,
  SNES_NCG_PRP = 1,
  SNES_NCG_HS = 2,
  SNES_NCG_DY = 3,
  SNES_NCG_CD = 4} SNESNCGType;
extern const char *const SNESNCGTypes[];

extern PetscErrorCode SNESNCGSetType(SNES, SNESNCGType);

typedef enum {SNES_QN_SCALE_DEFAULT = 0,
              SNES_QN_SCALE_NONE = 1,
              SNES_QN_SCALE_SHANNO = 2,
              SNES_QN_SCALE_LINESEARCH = 3,
              SNES_QN_SCALE_JACOBIAN = 4} SNESQNScaleType;
extern const char *const SNESQNScaleTypes[];
typedef enum {SNES_QN_RESTART_DEFAULT = 0,
              SNES_QN_RESTART_NONE = 1,
              SNES_QN_RESTART_POWELL = 2,
              SNES_QN_RESTART_PERIODIC = 3} SNESQNRestartType;
extern const char *const SNESQNRestartTypes[];
typedef enum {SNES_QN_LBFGS = 0,
              SNES_QN_BROYDEN = 1,
              SNES_QN_BADBROYDEN = 2
             } SNESQNType;
extern const char *const SNESQNTypes[];

extern PetscErrorCode SNESQNSetType(SNES, SNESQNType);
extern PetscErrorCode SNESQNSetScaleType(SNES, SNESQNScaleType);
extern PetscErrorCode SNESQNSetRestartType(SNES, SNESQNRestartType);

extern PetscErrorCode SNESNASMGetType(SNES,PCASMType*);
extern PetscErrorCode SNESNASMSetType(SNES,PCASMType);
extern PetscErrorCode SNESNASMGetSubdomains(SNES,PetscInt*,SNES**,VecScatter**,VecScatter**,VecScatter**);
extern PetscErrorCode SNESNASMSetSubdomains(SNES,PetscInt,SNES*,VecScatter*,VecScatter*,VecScatter*);
extern PetscErrorCode SNESNASMSetDamping(SNES,PetscReal);
extern PetscErrorCode SNESNASMGetDamping(SNES,PetscReal*);
extern PetscErrorCode SNESNASMGetSubdomainVecs(SNES,PetscInt*,Vec**,Vec**,Vec**,Vec**);
extern PetscErrorCode SNESNASMSetComputeFinalJacobian(SNES,PetscBool);

typedef enum {SNES_COMPOSITE_ADDITIVE,SNES_COMPOSITE_MULTIPLICATIVE,SNES_COMPOSITE_ADDITIVEOPTIMAL} SNESCompositeType;
extern const char *const SNESCompositeTypes[];

extern PetscErrorCode SNESCompositeSetType(SNES,SNESCompositeType);
extern PetscErrorCode SNESCompositeAddSNES(SNES,SNESType);
extern PetscErrorCode SNESCompositeGetSNES(SNES,PetscInt,SNES *);
extern PetscErrorCode SNESCompositeGetNumber(SNES,PetscInt*);
extern PetscErrorCode SNESCompositeSetDamping(SNES,PetscInt,PetscReal);
# 785 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscsnes.h"
typedef enum { SNES_FAS_MULTIPLICATIVE, SNES_FAS_ADDITIVE, SNES_FAS_FULL, SNES_FAS_KASKADE } SNESFASType;
extern const char *const SNESFASTypes[];


extern PetscErrorCode SNESFASSetType(SNES, SNESFASType);
extern PetscErrorCode SNESFASGetType(SNES, SNESFASType*);
extern PetscErrorCode SNESFASSetLevels(SNES, PetscInt, MPI_Comm *);
extern PetscErrorCode SNESFASGetLevels(SNES, PetscInt *);
extern PetscErrorCode SNESFASGetCycleSNES(SNES, PetscInt, SNES*);
extern PetscErrorCode SNESFASSetNumberSmoothUp(SNES, PetscInt);
extern PetscErrorCode SNESFASSetNumberSmoothDown(SNES, PetscInt);
extern PetscErrorCode SNESFASSetCycles(SNES, PetscInt);
extern PetscErrorCode SNESFASSetMonitor(SNES, PetscBool);
extern PetscErrorCode SNESFASSetLog(SNES, PetscBool);

extern PetscErrorCode SNESFASSetGalerkin(SNES, PetscBool);
extern PetscErrorCode SNESFASGetGalerkin(SNES, PetscBool*);


extern PetscErrorCode SNESFASCycleGetSmoother(SNES, SNES*);
extern PetscErrorCode SNESFASCycleGetSmootherUp(SNES, SNES*);
extern PetscErrorCode SNESFASCycleGetSmootherDown(SNES, SNES*);
extern PetscErrorCode SNESFASCycleGetCorrection(SNES, SNES*);
extern PetscErrorCode SNESFASCycleGetInterpolation(SNES, Mat*);
extern PetscErrorCode SNESFASCycleGetRestriction(SNES, Mat*);
extern PetscErrorCode SNESFASCycleGetInjection(SNES, Mat*);
extern PetscErrorCode SNESFASCycleGetRScale(SNES, Vec*);
extern PetscErrorCode SNESFASCycleSetCycles(SNES, PetscInt);
extern PetscErrorCode SNESFASCycleIsFine(SNES, PetscBool*);


extern PetscErrorCode SNESFASSetInterpolation(SNES, PetscInt, Mat);
extern PetscErrorCode SNESFASGetInterpolation(SNES, PetscInt, Mat*);
extern PetscErrorCode SNESFASSetRestriction(SNES, PetscInt, Mat);
extern PetscErrorCode SNESFASGetRestriction(SNES, PetscInt, Mat*);
extern PetscErrorCode SNESFASSetInjection(SNES, PetscInt, Mat);
extern PetscErrorCode SNESFASGetInjection(SNES, PetscInt, Mat*);
extern PetscErrorCode SNESFASSetRScale(SNES, PetscInt, Vec);
extern PetscErrorCode SNESFASGetRScale(SNES, PetscInt, Vec*);
extern PetscErrorCode SNESFASSetContinuation(SNES,PetscBool);

extern PetscErrorCode SNESFASGetSmoother(SNES, PetscInt, SNES*);
extern PetscErrorCode SNESFASGetSmootherUp(SNES, PetscInt, SNES*);
extern PetscErrorCode SNESFASGetSmootherDown(SNES, PetscInt, SNES*);
extern PetscErrorCode SNESFASGetCoarseSolve(SNES, SNES*);


extern PetscErrorCode SNESFASFullSetDownSweep(SNES,PetscBool);
extern PetscErrorCode SNESFASCreateCoarseVec(SNES,Vec*);
extern PetscErrorCode SNESFASRestrict(SNES,Vec,Vec);

extern PetscErrorCode DMSNESCheckFromOptions(SNES,Vec,PetscErrorCode (**)(PetscInt,const PetscReal[],PetscInt,PetscScalar*,void*),void**);
# 5 "variables.h" 2
# 1 "/usr/include/unistd.h" 1 3 4
# 27 "/usr/include/unistd.h" 3 4

# 202 "/usr/include/unistd.h" 3 4
# 1 "/usr/include/bits/posix_opt.h" 1 3 4
# 203 "/usr/include/unistd.h" 2 3 4



# 1 "/usr/include/bits/environments.h" 1 3 4
# 22 "/usr/include/bits/environments.h" 3 4
# 1 "/usr/include/bits/wordsize.h" 1 3 4
# 23 "/usr/include/bits/environments.h" 2 3 4
# 207 "/usr/include/unistd.h" 2 3 4
# 226 "/usr/include/unistd.h" 3 4
# 1 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/stddef.h" 1 3 4
# 227 "/usr/include/unistd.h" 2 3 4
# 255 "/usr/include/unistd.h" 3 4

# 255 "/usr/include/unistd.h" 3 4
typedef __useconds_t useconds_t;
# 274 "/usr/include/unistd.h" 3 4
typedef __socklen_t socklen_t;
# 287 "/usr/include/unistd.h" 3 4
extern int access (const char *__name, int __type) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));
# 304 "/usr/include/unistd.h" 3 4
extern int faccessat (int __fd, const char *__file, int __type, int __flag)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (2))) ;
# 334 "/usr/include/unistd.h" 3 4
extern __off_t lseek (int __fd, __off_t __offset, int __whence) __attribute__ ((__nothrow__ , __leaf__));
# 353 "/usr/include/unistd.h" 3 4
extern int close (int __fd);






extern ssize_t read (int __fd, void *__buf, size_t __nbytes) ;





extern ssize_t write (int __fd, const void *__buf, size_t __n) ;
# 376 "/usr/include/unistd.h" 3 4
extern ssize_t pread (int __fd, void *__buf, size_t __nbytes,
        __off_t __offset) ;






extern ssize_t pwrite (int __fd, const void *__buf, size_t __n,
         __off_t __offset) ;
# 417 "/usr/include/unistd.h" 3 4
extern int pipe (int __pipedes[2]) __attribute__ ((__nothrow__ , __leaf__)) ;
# 432 "/usr/include/unistd.h" 3 4
extern unsigned int alarm (unsigned int __seconds) __attribute__ ((__nothrow__ , __leaf__));
# 444 "/usr/include/unistd.h" 3 4
extern unsigned int sleep (unsigned int __seconds);







extern __useconds_t ualarm (__useconds_t __value, __useconds_t __interval)
     __attribute__ ((__nothrow__ , __leaf__));






extern int usleep (__useconds_t __useconds);
# 469 "/usr/include/unistd.h" 3 4
extern int pause (void);



extern int chown (const char *__file, __uid_t __owner, __gid_t __group)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1))) ;



extern int fchown (int __fd, __uid_t __owner, __gid_t __group) __attribute__ ((__nothrow__ , __leaf__)) ;




extern int lchown (const char *__file, __uid_t __owner, __gid_t __group)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1))) ;






extern int fchownat (int __fd, const char *__file, __uid_t __owner,
       __gid_t __group, int __flag)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (2))) ;



extern int chdir (const char *__path) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1))) ;



extern int fchdir (int __fd) __attribute__ ((__nothrow__ , __leaf__)) ;
# 511 "/usr/include/unistd.h" 3 4
extern char *getcwd (char *__buf, size_t __size) __attribute__ ((__nothrow__ , __leaf__)) ;
# 525 "/usr/include/unistd.h" 3 4
extern char *getwd (char *__buf)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1))) __attribute__ ((__deprecated__)) ;




extern int dup (int __fd) __attribute__ ((__nothrow__ , __leaf__)) ;


extern int dup2 (int __fd, int __fd2) __attribute__ ((__nothrow__ , __leaf__));
# 543 "/usr/include/unistd.h" 3 4
extern char **__environ;







extern int execve (const char *__path, char *const __argv[],
     char *const __envp[]) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));




extern int fexecve (int __fd, char *const __argv[], char *const __envp[])
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (2)));




extern int execv (const char *__path, char *const __argv[])
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));



extern int execle (const char *__path, const char *__arg, ...)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));



extern int execl (const char *__path, const char *__arg, ...)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));



extern int execvp (const char *__file, char *const __argv[])
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));




extern int execlp (const char *__file, const char *__arg, ...)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));
# 598 "/usr/include/unistd.h" 3 4
extern int nice (int __inc) __attribute__ ((__nothrow__ , __leaf__)) ;




extern void _exit (int __status) __attribute__ ((__noreturn__));





# 1 "/usr/include/bits/confname.h" 1 3 4
# 25 "/usr/include/bits/confname.h" 3 4
enum
  {
    _PC_LINK_MAX,

    _PC_MAX_CANON,

    _PC_MAX_INPUT,

    _PC_NAME_MAX,

    _PC_PATH_MAX,

    _PC_PIPE_BUF,

    _PC_CHOWN_RESTRICTED,

    _PC_NO_TRUNC,

    _PC_VDISABLE,

    _PC_SYNC_IO,

    _PC_ASYNC_IO,

    _PC_PRIO_IO,

    _PC_SOCK_MAXBUF,

    _PC_FILESIZEBITS,

    _PC_REC_INCR_XFER_SIZE,

    _PC_REC_MAX_XFER_SIZE,

    _PC_REC_MIN_XFER_SIZE,

    _PC_REC_XFER_ALIGN,

    _PC_ALLOC_SIZE_MIN,

    _PC_SYMLINK_MAX,

    _PC_2_SYMLINKS

  };


enum
  {
    _SC_ARG_MAX,

    _SC_CHILD_MAX,

    _SC_CLK_TCK,

    _SC_NGROUPS_MAX,

    _SC_OPEN_MAX,

    _SC_STREAM_MAX,

    _SC_TZNAME_MAX,

    _SC_JOB_CONTROL,

    _SC_SAVED_IDS,

    _SC_REALTIME_SIGNALS,

    _SC_PRIORITY_SCHEDULING,

    _SC_TIMERS,

    _SC_ASYNCHRONOUS_IO,

    _SC_PRIORITIZED_IO,

    _SC_SYNCHRONIZED_IO,

    _SC_FSYNC,

    _SC_MAPPED_FILES,

    _SC_MEMLOCK,

    _SC_MEMLOCK_RANGE,

    _SC_MEMORY_PROTECTION,

    _SC_MESSAGE_PASSING,

    _SC_SEMAPHORES,

    _SC_SHARED_MEMORY_OBJECTS,

    _SC_AIO_LISTIO_MAX,

    _SC_AIO_MAX,

    _SC_AIO_PRIO_DELTA_MAX,

    _SC_DELAYTIMER_MAX,

    _SC_MQ_OPEN_MAX,

    _SC_MQ_PRIO_MAX,

    _SC_VERSION,

    _SC_PAGESIZE,


    _SC_RTSIG_MAX,

    _SC_SEM_NSEMS_MAX,

    _SC_SEM_VALUE_MAX,

    _SC_SIGQUEUE_MAX,

    _SC_TIMER_MAX,




    _SC_BC_BASE_MAX,

    _SC_BC_DIM_MAX,

    _SC_BC_SCALE_MAX,

    _SC_BC_STRING_MAX,

    _SC_COLL_WEIGHTS_MAX,

    _SC_EQUIV_CLASS_MAX,

    _SC_EXPR_NEST_MAX,

    _SC_LINE_MAX,

    _SC_RE_DUP_MAX,

    _SC_CHARCLASS_NAME_MAX,


    _SC_2_VERSION,

    _SC_2_C_BIND,

    _SC_2_C_DEV,

    _SC_2_FORT_DEV,

    _SC_2_FORT_RUN,

    _SC_2_SW_DEV,

    _SC_2_LOCALEDEF,


    _SC_PII,

    _SC_PII_XTI,

    _SC_PII_SOCKET,

    _SC_PII_INTERNET,

    _SC_PII_OSI,

    _SC_POLL,

    _SC_SELECT,

    _SC_UIO_MAXIOV,

    _SC_IOV_MAX = _SC_UIO_MAXIOV,

    _SC_PII_INTERNET_STREAM,

    _SC_PII_INTERNET_DGRAM,

    _SC_PII_OSI_COTS,

    _SC_PII_OSI_CLTS,

    _SC_PII_OSI_M,

    _SC_T_IOV_MAX,



    _SC_THREADS,

    _SC_THREAD_SAFE_FUNCTIONS,

    _SC_GETGR_R_SIZE_MAX,

    _SC_GETPW_R_SIZE_MAX,

    _SC_LOGIN_NAME_MAX,

    _SC_TTY_NAME_MAX,

    _SC_THREAD_DESTRUCTOR_ITERATIONS,

    _SC_THREAD_KEYS_MAX,

    _SC_THREAD_STACK_MIN,

    _SC_THREAD_THREADS_MAX,

    _SC_THREAD_ATTR_STACKADDR,

    _SC_THREAD_ATTR_STACKSIZE,

    _SC_THREAD_PRIORITY_SCHEDULING,

    _SC_THREAD_PRIO_INHERIT,

    _SC_THREAD_PRIO_PROTECT,

    _SC_THREAD_PROCESS_SHARED,


    _SC_NPROCESSORS_CONF,

    _SC_NPROCESSORS_ONLN,

    _SC_PHYS_PAGES,

    _SC_AVPHYS_PAGES,

    _SC_ATEXIT_MAX,

    _SC_PASS_MAX,


    _SC_XOPEN_VERSION,

    _SC_XOPEN_XCU_VERSION,

    _SC_XOPEN_UNIX,

    _SC_XOPEN_CRYPT,

    _SC_XOPEN_ENH_I18N,

    _SC_XOPEN_SHM,


    _SC_2_CHAR_TERM,

    _SC_2_C_VERSION,

    _SC_2_UPE,


    _SC_XOPEN_XPG2,

    _SC_XOPEN_XPG3,

    _SC_XOPEN_XPG4,


    _SC_CHAR_BIT,

    _SC_CHAR_MAX,

    _SC_CHAR_MIN,

    _SC_INT_MAX,

    _SC_INT_MIN,

    _SC_LONG_BIT,

    _SC_WORD_BIT,

    _SC_MB_LEN_MAX,

    _SC_NZERO,

    _SC_SSIZE_MAX,

    _SC_SCHAR_MAX,

    _SC_SCHAR_MIN,

    _SC_SHRT_MAX,

    _SC_SHRT_MIN,

    _SC_UCHAR_MAX,

    _SC_UINT_MAX,

    _SC_ULONG_MAX,

    _SC_USHRT_MAX,


    _SC_NL_ARGMAX,

    _SC_NL_LANGMAX,

    _SC_NL_MSGMAX,

    _SC_NL_NMAX,

    _SC_NL_SETMAX,

    _SC_NL_TEXTMAX,


    _SC_XBS5_ILP32_OFF32,

    _SC_XBS5_ILP32_OFFBIG,

    _SC_XBS5_LP64_OFF64,

    _SC_XBS5_LPBIG_OFFBIG,


    _SC_XOPEN_LEGACY,

    _SC_XOPEN_REALTIME,

    _SC_XOPEN_REALTIME_THREADS,


    _SC_ADVISORY_INFO,

    _SC_BARRIERS,

    _SC_BASE,

    _SC_C_LANG_SUPPORT,

    _SC_C_LANG_SUPPORT_R,

    _SC_CLOCK_SELECTION,

    _SC_CPUTIME,

    _SC_THREAD_CPUTIME,

    _SC_DEVICE_IO,

    _SC_DEVICE_SPECIFIC,

    _SC_DEVICE_SPECIFIC_R,

    _SC_FD_MGMT,

    _SC_FIFO,

    _SC_PIPE,

    _SC_FILE_ATTRIBUTES,

    _SC_FILE_LOCKING,

    _SC_FILE_SYSTEM,

    _SC_MONOTONIC_CLOCK,

    _SC_MULTI_PROCESS,

    _SC_SINGLE_PROCESS,

    _SC_NETWORKING,

    _SC_READER_WRITER_LOCKS,

    _SC_SPIN_LOCKS,

    _SC_REGEXP,

    _SC_REGEX_VERSION,

    _SC_SHELL,

    _SC_SIGNALS,

    _SC_SPAWN,

    _SC_SPORADIC_SERVER,

    _SC_THREAD_SPORADIC_SERVER,

    _SC_SYSTEM_DATABASE,

    _SC_SYSTEM_DATABASE_R,

    _SC_TIMEOUTS,

    _SC_TYPED_MEMORY_OBJECTS,

    _SC_USER_GROUPS,

    _SC_USER_GROUPS_R,

    _SC_2_PBS,

    _SC_2_PBS_ACCOUNTING,

    _SC_2_PBS_LOCATE,

    _SC_2_PBS_MESSAGE,

    _SC_2_PBS_TRACK,

    _SC_SYMLOOP_MAX,

    _SC_STREAMS,

    _SC_2_PBS_CHECKPOINT,


    _SC_V6_ILP32_OFF32,

    _SC_V6_ILP32_OFFBIG,

    _SC_V6_LP64_OFF64,

    _SC_V6_LPBIG_OFFBIG,


    _SC_HOST_NAME_MAX,

    _SC_TRACE,

    _SC_TRACE_EVENT_FILTER,

    _SC_TRACE_INHERIT,

    _SC_TRACE_LOG,


    _SC_LEVEL1_ICACHE_SIZE,

    _SC_LEVEL1_ICACHE_ASSOC,

    _SC_LEVEL1_ICACHE_LINESIZE,

    _SC_LEVEL1_DCACHE_SIZE,

    _SC_LEVEL1_DCACHE_ASSOC,

    _SC_LEVEL1_DCACHE_LINESIZE,

    _SC_LEVEL2_CACHE_SIZE,

    _SC_LEVEL2_CACHE_ASSOC,

    _SC_LEVEL2_CACHE_LINESIZE,

    _SC_LEVEL3_CACHE_SIZE,

    _SC_LEVEL3_CACHE_ASSOC,

    _SC_LEVEL3_CACHE_LINESIZE,

    _SC_LEVEL4_CACHE_SIZE,

    _SC_LEVEL4_CACHE_ASSOC,

    _SC_LEVEL4_CACHE_LINESIZE,



    _SC_IPV6 = _SC_LEVEL1_ICACHE_SIZE + 50,

    _SC_RAW_SOCKETS,


    _SC_V7_ILP32_OFF32,

    _SC_V7_ILP32_OFFBIG,

    _SC_V7_LP64_OFF64,

    _SC_V7_LPBIG_OFFBIG,


    _SC_SS_REPL_MAX,


    _SC_TRACE_EVENT_NAME_MAX,

    _SC_TRACE_NAME_MAX,

    _SC_TRACE_SYS_MAX,

    _SC_TRACE_USER_EVENT_MAX,


    _SC_XOPEN_STREAMS,


    _SC_THREAD_ROBUST_PRIO_INHERIT,

    _SC_THREAD_ROBUST_PRIO_PROTECT

  };


enum
  {
    _CS_PATH,


    _CS_V6_WIDTH_RESTRICTED_ENVS,



    _CS_GNU_LIBC_VERSION,

    _CS_GNU_LIBPTHREAD_VERSION,


    _CS_V5_WIDTH_RESTRICTED_ENVS,



    _CS_V7_WIDTH_RESTRICTED_ENVS,



    _CS_LFS_CFLAGS = 1000,

    _CS_LFS_LDFLAGS,

    _CS_LFS_LIBS,

    _CS_LFS_LINTFLAGS,

    _CS_LFS64_CFLAGS,

    _CS_LFS64_LDFLAGS,

    _CS_LFS64_LIBS,

    _CS_LFS64_LINTFLAGS,


    _CS_XBS5_ILP32_OFF32_CFLAGS = 1100,

    _CS_XBS5_ILP32_OFF32_LDFLAGS,

    _CS_XBS5_ILP32_OFF32_LIBS,

    _CS_XBS5_ILP32_OFF32_LINTFLAGS,

    _CS_XBS5_ILP32_OFFBIG_CFLAGS,

    _CS_XBS5_ILP32_OFFBIG_LDFLAGS,

    _CS_XBS5_ILP32_OFFBIG_LIBS,

    _CS_XBS5_ILP32_OFFBIG_LINTFLAGS,

    _CS_XBS5_LP64_OFF64_CFLAGS,

    _CS_XBS5_LP64_OFF64_LDFLAGS,

    _CS_XBS5_LP64_OFF64_LIBS,

    _CS_XBS5_LP64_OFF64_LINTFLAGS,

    _CS_XBS5_LPBIG_OFFBIG_CFLAGS,

    _CS_XBS5_LPBIG_OFFBIG_LDFLAGS,

    _CS_XBS5_LPBIG_OFFBIG_LIBS,

    _CS_XBS5_LPBIG_OFFBIG_LINTFLAGS,


    _CS_POSIX_V6_ILP32_OFF32_CFLAGS,

    _CS_POSIX_V6_ILP32_OFF32_LDFLAGS,

    _CS_POSIX_V6_ILP32_OFF32_LIBS,

    _CS_POSIX_V6_ILP32_OFF32_LINTFLAGS,

    _CS_POSIX_V6_ILP32_OFFBIG_CFLAGS,

    _CS_POSIX_V6_ILP32_OFFBIG_LDFLAGS,

    _CS_POSIX_V6_ILP32_OFFBIG_LIBS,

    _CS_POSIX_V6_ILP32_OFFBIG_LINTFLAGS,

    _CS_POSIX_V6_LP64_OFF64_CFLAGS,

    _CS_POSIX_V6_LP64_OFF64_LDFLAGS,

    _CS_POSIX_V6_LP64_OFF64_LIBS,

    _CS_POSIX_V6_LP64_OFF64_LINTFLAGS,

    _CS_POSIX_V6_LPBIG_OFFBIG_CFLAGS,

    _CS_POSIX_V6_LPBIG_OFFBIG_LDFLAGS,

    _CS_POSIX_V6_LPBIG_OFFBIG_LIBS,

    _CS_POSIX_V6_LPBIG_OFFBIG_LINTFLAGS,


    _CS_POSIX_V7_ILP32_OFF32_CFLAGS,

    _CS_POSIX_V7_ILP32_OFF32_LDFLAGS,

    _CS_POSIX_V7_ILP32_OFF32_LIBS,

    _CS_POSIX_V7_ILP32_OFF32_LINTFLAGS,

    _CS_POSIX_V7_ILP32_OFFBIG_CFLAGS,

    _CS_POSIX_V7_ILP32_OFFBIG_LDFLAGS,

    _CS_POSIX_V7_ILP32_OFFBIG_LIBS,

    _CS_POSIX_V7_ILP32_OFFBIG_LINTFLAGS,

    _CS_POSIX_V7_LP64_OFF64_CFLAGS,

    _CS_POSIX_V7_LP64_OFF64_LDFLAGS,

    _CS_POSIX_V7_LP64_OFF64_LIBS,

    _CS_POSIX_V7_LP64_OFF64_LINTFLAGS,

    _CS_POSIX_V7_LPBIG_OFFBIG_CFLAGS,

    _CS_POSIX_V7_LPBIG_OFFBIG_LDFLAGS,

    _CS_POSIX_V7_LPBIG_OFFBIG_LIBS,

    _CS_POSIX_V7_LPBIG_OFFBIG_LINTFLAGS,


    _CS_V6_ENV,

    _CS_V7_ENV

  };
# 610 "/usr/include/unistd.h" 2 3 4


extern long int pathconf (const char *__path, int __name)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));


extern long int fpathconf (int __fd, int __name) __attribute__ ((__nothrow__ , __leaf__));


extern long int sysconf (int __name) __attribute__ ((__nothrow__ , __leaf__));



extern size_t confstr (int __name, char *__buf, size_t __len) __attribute__ ((__nothrow__ , __leaf__));




extern __pid_t getpid (void) __attribute__ ((__nothrow__ , __leaf__));


extern __pid_t getppid (void) __attribute__ ((__nothrow__ , __leaf__));




extern __pid_t getpgrp (void) __attribute__ ((__nothrow__ , __leaf__));
# 646 "/usr/include/unistd.h" 3 4
extern __pid_t __getpgid (__pid_t __pid) __attribute__ ((__nothrow__ , __leaf__));

extern __pid_t getpgid (__pid_t __pid) __attribute__ ((__nothrow__ , __leaf__));






extern int setpgid (__pid_t __pid, __pid_t __pgid) __attribute__ ((__nothrow__ , __leaf__));
# 672 "/usr/include/unistd.h" 3 4
extern int setpgrp (void) __attribute__ ((__nothrow__ , __leaf__));
# 689 "/usr/include/unistd.h" 3 4
extern __pid_t setsid (void) __attribute__ ((__nothrow__ , __leaf__));



extern __pid_t getsid (__pid_t __pid) __attribute__ ((__nothrow__ , __leaf__));



extern __uid_t getuid (void) __attribute__ ((__nothrow__ , __leaf__));


extern __uid_t geteuid (void) __attribute__ ((__nothrow__ , __leaf__));


extern __gid_t getgid (void) __attribute__ ((__nothrow__ , __leaf__));


extern __gid_t getegid (void) __attribute__ ((__nothrow__ , __leaf__));




extern int getgroups (int __size, __gid_t __list[]) __attribute__ ((__nothrow__ , __leaf__)) ;
# 722 "/usr/include/unistd.h" 3 4
extern int setuid (__uid_t __uid) __attribute__ ((__nothrow__ , __leaf__)) ;




extern int setreuid (__uid_t __ruid, __uid_t __euid) __attribute__ ((__nothrow__ , __leaf__)) ;




extern int seteuid (__uid_t __uid) __attribute__ ((__nothrow__ , __leaf__)) ;






extern int setgid (__gid_t __gid) __attribute__ ((__nothrow__ , __leaf__)) ;




extern int setregid (__gid_t __rgid, __gid_t __egid) __attribute__ ((__nothrow__ , __leaf__)) ;




extern int setegid (__gid_t __gid) __attribute__ ((__nothrow__ , __leaf__)) ;
# 778 "/usr/include/unistd.h" 3 4
extern __pid_t fork (void) __attribute__ ((__nothrow__));







extern __pid_t vfork (void) __attribute__ ((__nothrow__ , __leaf__));





extern char *ttyname (int __fd) __attribute__ ((__nothrow__ , __leaf__));



extern int ttyname_r (int __fd, char *__buf, size_t __buflen)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (2))) ;



extern int isatty (int __fd) __attribute__ ((__nothrow__ , __leaf__));





extern int ttyslot (void) __attribute__ ((__nothrow__ , __leaf__));




extern int link (const char *__from, const char *__to)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2))) ;




extern int linkat (int __fromfd, const char *__from, int __tofd,
     const char *__to, int __flags)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (2, 4))) ;




extern int symlink (const char *__from, const char *__to)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2))) ;




extern ssize_t readlink (const char *__restrict __path,
    char *__restrict __buf, size_t __len)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2))) ;




extern int symlinkat (const char *__from, int __tofd,
        const char *__to) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 3))) ;


extern ssize_t readlinkat (int __fd, const char *__restrict __path,
      char *__restrict __buf, size_t __len)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (2, 3))) ;



extern int unlink (const char *__name) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));



extern int unlinkat (int __fd, const char *__name, int __flag)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (2)));



extern int rmdir (const char *__path) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));



extern __pid_t tcgetpgrp (int __fd) __attribute__ ((__nothrow__ , __leaf__));


extern int tcsetpgrp (int __fd, __pid_t __pgrp_id) __attribute__ ((__nothrow__ , __leaf__));






extern char *getlogin (void);







extern int getlogin_r (char *__name, size_t __name_len) __attribute__ ((__nonnull__ (1)));




extern int setlogin (const char *__name) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));
# 893 "/usr/include/unistd.h" 3 4
# 1 "/usr/include/getopt.h" 1 3 4
# 58 "/usr/include/getopt.h" 3 4
extern char *optarg;
# 72 "/usr/include/getopt.h" 3 4
extern int optind;




extern int opterr;



extern int optopt;
# 151 "/usr/include/getopt.h" 3 4
extern int getopt (int ___argc, char *const *___argv, const char *__shortopts)
       __attribute__ ((__nothrow__ , __leaf__));
# 894 "/usr/include/unistd.h" 2 3 4







extern int gethostname (char *__name, size_t __len) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));






extern int sethostname (const char *__name, size_t __len)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1))) ;



extern int sethostid (long int __id) __attribute__ ((__nothrow__ , __leaf__)) ;





extern int getdomainname (char *__name, size_t __len)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1))) ;
extern int setdomainname (const char *__name, size_t __len)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1))) ;





extern int vhangup (void) __attribute__ ((__nothrow__ , __leaf__));


extern int revoke (const char *__file) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1))) ;







extern int profil (unsigned short int *__sample_buffer, size_t __size,
     size_t __offset, unsigned int __scale)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));





extern int acct (const char *__name) __attribute__ ((__nothrow__ , __leaf__));



extern char *getusershell (void) __attribute__ ((__nothrow__ , __leaf__));
extern void endusershell (void) __attribute__ ((__nothrow__ , __leaf__));
extern void setusershell (void) __attribute__ ((__nothrow__ , __leaf__));





extern int daemon (int __nochdir, int __noclose) __attribute__ ((__nothrow__ , __leaf__)) ;






extern int chroot (const char *__path) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1))) ;



extern char *getpass (const char *__prompt) __attribute__ ((__nonnull__ (1)));







extern int fsync (int __fd);
# 991 "/usr/include/unistd.h" 3 4
extern long int gethostid (void);


extern void sync (void) __attribute__ ((__nothrow__ , __leaf__));





extern int getpagesize (void) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));




extern int getdtablesize (void) __attribute__ ((__nothrow__ , __leaf__));
# 1015 "/usr/include/unistd.h" 3 4
extern int truncate (const char *__file, __off_t __length)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1))) ;
# 1038 "/usr/include/unistd.h" 3 4
extern int ftruncate (int __fd, __off_t __length) __attribute__ ((__nothrow__ , __leaf__)) ;
# 1059 "/usr/include/unistd.h" 3 4
extern int brk (void *__addr) __attribute__ ((__nothrow__ , __leaf__)) ;





extern void *sbrk (intptr_t __delta) __attribute__ ((__nothrow__ , __leaf__));
# 1080 "/usr/include/unistd.h" 3 4
extern long int syscall (long int __sysno, ...) __attribute__ ((__nothrow__ , __leaf__));
# 1103 "/usr/include/unistd.h" 3 4
extern int lockf (int __fd, int __cmd, __off_t __len) ;
# 1134 "/usr/include/unistd.h" 3 4
extern int fdatasync (int __fildes);
# 1172 "/usr/include/unistd.h" 3 4

# 6 "variables.h" 2
# 18 "variables.h"

# 18 "variables.h"
typedef struct {
  PetscScalar x, y, z;
} Cmpnts;

typedef struct {
  PetscReal t, f;
} FlowWave;

typedef struct {
 PetscScalar x, y;
} Cmpnts2;

typedef struct {
  PetscReal x, y;
} Cpt2D;

typedef struct {
  Vec Ubcs;
  } BCS;

typedef struct {
  Vec Res, x, Fext, Fint, Fdyn, disp, FJ, dis, V;
  Vec xn, xnm1, xd, xdsc, xdd, xddsc, dx, dxn, dxsc, xnsc, y, yn, xninner, yninner;

  PetscReal *StrainM, *StressM, *StrainB, *StressB, *Pf, *Pfn;
  PetscReal dt, wAit;

  Vec Mass, Dissip;

    Vec Fcnt;
    Cmpnts *qvec;
    PetscReal *radvec;

    PetscInt *ire, *irv, *val, *patch;
    PetscReal *G, *G1, *G2;
    PetscInt *contact;

    PetscInt n_v, n_elmt, ibi, n_edge, sum_n_bnodes, n_ghosts;
    PetscReal *x_bp, *y_bp, *z_bp, *x_bp0, *y_bp0, *z_bp0, *x_bpr, *y_bpr, *z_bpr;
    PetscReal *p4x, *p4y, *p4z, *p5x, *p5y, *p5z, *p6x ,*p6y, *p6z;
    PetscReal *p4x0, *p4y0, *p4z0, *p5x0, *p5y0, *p5z0, *p6x0, *p6y0, *p6z0;
    PetscReal *kve0, *kve;
    PetscInt *nv1, *nv2, *nv3, *nv4, *nv5, *nv6;
    Cmpnts *n_fib;
    PetscInt *bnodes, *n_bnodes, *belmts, *edgefrontnodes, *edgefrontnodesI, *belmtsedge;
    PetscReal *nf_x, *nf_y, *nf_z, *Nf_x, *Nf_y, *Nf_z;
    PetscReal *dA, *dA0;
  } FE;

typedef struct {
  PetscInt i1, j1, k1;
  PetscInt i2, j2, k2;
  PetscInt i3, j3, k3;
  PetscReal cr1, cr2, cr3;
  PetscReal d_i;
  PetscInt imode;

  PetscInt ni, nj, nk;
  PetscReal d_s;
  Cmpnts pmin;
  PetscInt cell;
  PetscReal cs1, cs2, cs3;

  PetscInt i11, j11, k11;
  PetscInt i22, j22, k22;
  PetscInt i33, j33, k33;
  PetscReal cr11, cr22, cr33;
  PetscReal d_ii;
  PetscInt iimode;
  PetscReal cs11, cs22, cs33;

  PetscInt ii1, jj1, kk1;
  PetscInt ii2, jj2, kk2;
  PetscInt ii3, jj3, kk3;
  PetscReal ct1, ct2, ct3;

  PetscInt smode;

  PetscInt ii11, jj11, kk11;
  PetscInt ii22, jj22, kk22;
  PetscInt ii33, jj33, kk33;
  PetscReal ct11, ct22, ct33;
  PetscReal d_ss;
  PetscInt ssmode;
# 110 "variables.h"
} IBMInfo;

typedef struct {
  PetscInt nbnumber;
  PetscInt n_v, n_elmt;
  PetscInt *nv1, *nv2, *nv3;
  PetscReal *nf_x, *nf_y, *nf_z;
  PetscReal *x_bp, *y_bp, *z_bp;
  PetscReal *x_bp0, *y_bp0, *z_bp0;
  PetscReal *x_bp_o, *y_bp_o, *z_bp_o;
  PetscReal x_bp_in[101][3270], y_bp_in[101][3270], z_bp_in[101][3270];
  Cmpnts *u, *uold, *urm1;


  PetscReal *dA ;

  PetscReal *nt_x, *nt_y, *nt_z;
  PetscReal *ns_x, *ns_y, *ns_z;


  PetscReal *cent_x,*cent_y,*cent_z;


  PetscReal *pres, *tau0, *tauN;
  PetscReal *Bvel_u, *Bvel_v, *Bvel_w;
  PetscReal x_min,x_max,y_min,y_max,z_min,z_max;

  Cmpnts *qvec;
  PetscReal *radvec;
} IBMNodes;
typedef struct {
  PetscInt nbnumber;
  PetscInt n_v, n_elmt;
  PetscInt *nv1, *nv2, *nv3, *nv4;

  PetscReal *x_bp, *y_bp, *z_bp;
  PetscReal *x_bp0, *y_bp0, *z_bp0;
  PetscReal *x_bp_o, *y_bp_o, *z_bp_o;

  Cmpnts *u, *uold, *urm1;
  PetscReal V;

  PetscReal *dV0 ;

  PetscReal *cent_x,*cent_y,*cent_z;
  PetscReal x_c,y_c,z_c;
  PetscReal J[3][3];
  PetscReal I_inv[3][3];
} IBMVNodes;
typedef struct node{
  PetscInt Node;
  struct node *next;
} node;

typedef struct list{
  node *head;
} List;


typedef struct list_node {
  PetscInt index;
  struct list_node *next;
} Node_List;

typedef struct IBMListNode {
  IBMInfo ibm_intp;
  struct IBMListNode* next;
} IBMListNode;

typedef struct IBMList {
  IBMListNode *head;
} IBMList;


typedef struct UserCtx {
  DM da;



  DM fda;

  DM pda;
  DMDALocalInfo info;

  Vec Cent;
  Vec Centx,Centy,Centz;
  Vec Csi, Eta, Zet, Aj;
  Vec ICsi, IEta, IZet, IAj;
  Vec JCsi, JEta, JZet, JAj;
  Vec KCsi, KEta, KZet, KAj;

  Vec Ucont,Vcont;
  Vec Ucat, Wcat;
  Vec Ucat_o;
  Vec Ucont_o, Ucont_rm1, Rhs, dUcont, pUcont;
  Vec P, P_o;

  Vec Qnew;
  Vec Ql;
  Vec shrr;
  Vec ItfcQ;
  Vec nhostQ;
  Vec lItfcQ;
  PetscReal LVOUT;


  PetscReal cdisx,cdisy,cdisz;

  Vec Phi;
  Vec GridSpace;
  Vec Nvert;
  Vec Nvert_o;
  Vec Itfc;
  BCS Bcs;

  Vec lUcont, lUcat, lP, lPhi, lVcont,lWcat;
  Vec lUcont_o, lUcont_rm1;
  Vec lCsi, lEta, lZet, lAj;
  Vec lICsi, lIEta, lIZet, lIAj;
  Vec lJCsi, lJEta, lJZet, lJAj;
  Vec lKCsi, lKEta, lKZet, lKAj;
  Vec lGridSpace;
  Vec lNvert, lNvert_o, lNFace;
  Vec lCent;
  Vec lItfc;

  Vec lMAreaCsi, lMAreaEta, lMAreaZet;

  Vec inletU;

  Vec nhostU;


  Vec DUold;

  Vec Forcing;
  Vec Ucont_MG;

  Vec Dt, Nu_t, CS;

  AO ao;

  PetscReal ren;
  PetscReal dt;
  PetscReal st;


  PetscReal FluxIntpSum;


  Vec psuedot;
  PetscReal cfl, vnn;

  PetscReal r[101], tin[101], uinr[101][1001];

  PetscInt ip[10],jp[10],kp[10],dispx[10],dispy[10],dispz[10],dispnn[10];
  PetscReal *itfchostx, *itfchosty, *itfchostz;
  PetscReal FluxInSum, FluxOutSum, FluxIntfcSum;
  PetscReal AreaOutSum, AreaIntfcSum;



  PetscErrorCode aotopetsc;
  PetscBool assignedA;

  PetscInt _this;
  PetscInt *idx_from;

  PetscInt bctype[6],inttype[7];
  PetscInt itfcptsnumber;
  PetscInt *itfcI, *itfcJ, *itfcK;
  PetscInt *itfchostI, *itfchostJ, *itfchostK, *itfchostB;
  PetscInt IM, JM, KM;
  PetscReal Max_X,Max_Y,Max_Z,Min_X,Min_Y,Min_Z;

  PetscInt ibmnumber;
  IBMInfo *ibm_intp;
  Mat A, C;
  KSP ksp;

  IBMNodes *ibm;

  DM *da_f, *da_c;
  struct UserCtx *user_f, *user_c;
  Vec *lNvert_c;

  Vec B;
  Vec Rhsp, X, R;

  Mat MR, MP;
  MatNullSpace nullsp;



  PetscInt *KSKE;
  PetscBool multinullspace;

  IBMList *ibmlist;

  PetscInt thislevel, mglevels;

  PetscInt isc, jsc, ksc;
  PetscInt cgrid;

  FlowWave *inflow;
  PetscInt number_flowwave;


  Vec Ucat_sum;
  Vec Ucat_cross_sum;
  Vec Ucat_square_sum;
  Vec P_sum;
  Vec lSx, lSy, lSz, lS;
  Vec lLM, lMM, lNM;

  Vec lNu_t, lF1;
  Vec lCs;
  Vec K_Omega, lK_Omega, K_Omega_o, lK_Omega_o;
  Vec Distance;
  Vec lItfcKO;

  PetscInt rhs_count;
  Vec Gid, Gidm;
  Vec Phi2, B2, Ucont2;
  PetscInt local_Phi2_size, p_global_begin;
  PetscInt reduced_p_size;

  DM fda2;
  Vec lUstar;
  Cmpnts2 **komega_plane;


  SNES snes;
  Mat J,PJ,J2;
  Vec RFC;
  ISColoring iscoloring;
  MatFDColoring matfdcoloring;



  Vec ParticleVec, lParticleVec;

} UserCtx;

typedef struct {
  UserCtx *user;
  PetscInt thislevel;
  DM packer;

} MGCtx;

typedef struct {
  PetscInt mglevels;
  PetscInt thislevel;

  PetscBool isc, jsc, ksc;
  MGCtx *mgctx;

  DM packer;
  SNES snespacker;
} UserMG;

typedef struct {
  PetscReal P;
  PetscInt n_P;
  PetscReal Tow_ws, Tow_wt;
  PetscReal Tow_wn;

  PetscInt Clsnbpt_i,Clsnbpt_j,Clsnbpt_k;
  PetscInt icell,jcell,kcell;
  PetscInt FoundAroundcell;
  PetscInt Need3rdPoint;

} SurfElmtInfo;

typedef struct {
  PetscReal S_new[6],S_old[6],S_real[6],S_realm1[6];
  PetscReal S_ang_n[6],S_ang_o[6],S_ang_r[6],S_ang_rm1[6];
  PetscReal red_vel, damp, mu_s;
  PetscReal F_x,F_y,F_z, A_tot;
  PetscReal F_x_old,F_y_old,F_z_old;
  PetscReal F_x_real,F_y_real,F_z_real;
  PetscReal M_x,M_y,M_z;
  PetscReal M_x_old,M_y_old,M_z_old;
  PetscReal M_x_real,M_y_real,M_z_real;
  PetscReal M_x_rm2,M_y_rm2,M_z_rm2;
  PetscReal M_x_rm3,M_y_rm3,M_z_rm3;
  PetscReal x_c,y_c,z_c;
  PetscReal a_c[3];
  PetscReal Mdpdn_x, Mdpdn_y,Mdpdn_z;
  PetscReal Mdpdn_x_old, Mdpdn_y_old,Mdpdn_z_old;
  PetscReal Power;

  PetscReal clone;
  PetscInt pbc[3];
  PetscReal I_inv[3][3];
  PetscReal L_n[3],L_o[3],L_r[3];
  PetscReal alpha[3];
  PetscReal acc[3];
  PetscReal R[3][3],q[4],q_r[4];

  PetscReal dS[6],dS_o[6],atk,atk_o;


  SurfElmtInfo *elmtinfo;
  IBMInfo *fsi_intp;


  PetscReal Max_xbc,Min_xbc;
  PetscReal Max_ybc,Min_ybc;
  PetscReal Max_zbc,Min_zbc;


  PetscInt CV_ys,CV_ye,CV_zs,CV_ze;
} FSInfo;

typedef struct {
  PetscInt n_time, n_midp, n_subit;
  PetscReal *x_midp, *y_midp, *z_midp;
  PetscReal *x_com, *y_com, *head_ang;
  PetscReal *s1,*s2,*s3;



  PetscReal *st1,*st2,*st3;




  Mat Mphi;
  PetscReal xmin,xmax,ymin,ymax,zmin,zmax;
} Cstart;

typedef struct{
  PetscInt cell[3];
  Cmpnts loc, weights,vel;
} Particle;

typedef struct{
   Cmpnts min_coords, max_coords;
} BoundingBox;



extern PetscInt i_periodic, j_periodic, k_periodic;
extern PetscInt ii_periodic, jj_periodic, kk_periodic;
extern PetscInt les, wallfunction, rans;
extern PetscInt i_homo_filter, j_homo_filter, k_homo_filter;
extern PetscInt testfilter_ik, testfilter_1d;
extern PetscInt clark;
# 5 "interpolation.c" 2

# 1 "/usr/include/time.h" 1 3 4
# 29 "/usr/include/time.h" 3 4








# 1 "/sw/eb/sw/GCCcore/8.3.0/lib/gcc/x86_64-pc-linux-gnu/8.3.0/include/stddef.h" 1 3 4
# 38 "/usr/include/time.h" 2 3 4



# 1 "/usr/include/bits/time.h" 1 3 4
# 42 "/usr/include/time.h" 2 3 4
# 131 "/usr/include/time.h" 3 4



# 133 "/usr/include/time.h" 3 4
struct tm
{
  int tm_sec;
  int tm_min;
  int tm_hour;
  int tm_mday;
  int tm_mon;
  int tm_year;
  int tm_wday;
  int tm_yday;
  int tm_isdst;


  long int tm_gmtoff;
  const char *tm_zone;




};








struct itimerspec
  {
    struct timespec it_interval;
    struct timespec it_value;
  };


struct sigevent;
# 186 "/usr/include/time.h" 3 4



extern clock_t clock (void) __attribute__ ((__nothrow__ , __leaf__));


extern time_t time (time_t *__timer) __attribute__ ((__nothrow__ , __leaf__));


extern double difftime (time_t __time1, time_t __time0)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));


extern time_t mktime (struct tm *__tp) __attribute__ ((__nothrow__ , __leaf__));





extern size_t strftime (char *__restrict __s, size_t __maxsize,
   const char *__restrict __format,
   const struct tm *__restrict __tp) __attribute__ ((__nothrow__ , __leaf__));

# 223 "/usr/include/time.h" 3 4
extern size_t strftime_l (char *__restrict __s, size_t __maxsize,
     const char *__restrict __format,
     const struct tm *__restrict __tp,
     __locale_t __loc) __attribute__ ((__nothrow__ , __leaf__));
# 236 "/usr/include/time.h" 3 4



extern struct tm *gmtime (const time_t *__timer) __attribute__ ((__nothrow__ , __leaf__));



extern struct tm *localtime (const time_t *__timer) __attribute__ ((__nothrow__ , __leaf__));





extern struct tm *gmtime_r (const time_t *__restrict __timer,
       struct tm *__restrict __tp) __attribute__ ((__nothrow__ , __leaf__));



extern struct tm *localtime_r (const time_t *__restrict __timer,
          struct tm *__restrict __tp) __attribute__ ((__nothrow__ , __leaf__));





extern char *asctime (const struct tm *__tp) __attribute__ ((__nothrow__ , __leaf__));


extern char *ctime (const time_t *__timer) __attribute__ ((__nothrow__ , __leaf__));







extern char *asctime_r (const struct tm *__restrict __tp,
   char *__restrict __buf) __attribute__ ((__nothrow__ , __leaf__));


extern char *ctime_r (const time_t *__restrict __timer,
        char *__restrict __buf) __attribute__ ((__nothrow__ , __leaf__));




extern char *__tzname[2];
extern int __daylight;
extern long int __timezone;




extern char *tzname[2];



extern void tzset (void) __attribute__ ((__nothrow__ , __leaf__));



extern int daylight;
extern long int timezone;





extern int stime (const time_t *__when) __attribute__ ((__nothrow__ , __leaf__));
# 319 "/usr/include/time.h" 3 4
extern time_t timegm (struct tm *__tp) __attribute__ ((__nothrow__ , __leaf__));


extern time_t timelocal (struct tm *__tp) __attribute__ ((__nothrow__ , __leaf__));


extern int dysize (int __year) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
# 334 "/usr/include/time.h" 3 4
extern int nanosleep (const struct timespec *__requested_time,
        struct timespec *__remaining);



extern int clock_getres (clockid_t __clock_id, struct timespec *__res) __attribute__ ((__nothrow__ , __leaf__));


extern int clock_gettime (clockid_t __clock_id, struct timespec *__tp) __attribute__ ((__nothrow__ , __leaf__));


extern int clock_settime (clockid_t __clock_id, const struct timespec *__tp)
     __attribute__ ((__nothrow__ , __leaf__));






extern int clock_nanosleep (clockid_t __clock_id, int __flags,
       const struct timespec *__req,
       struct timespec *__rem);


extern int clock_getcpuclockid (pid_t __pid, clockid_t *__clock_id) __attribute__ ((__nothrow__ , __leaf__));




extern int timer_create (clockid_t __clock_id,
    struct sigevent *__restrict __evp,
    timer_t *__restrict __timerid) __attribute__ ((__nothrow__ , __leaf__));


extern int timer_delete (timer_t __timerid) __attribute__ ((__nothrow__ , __leaf__));


extern int timer_settime (timer_t __timerid, int __flags,
     const struct itimerspec *__restrict __value,
     struct itimerspec *__restrict __ovalue) __attribute__ ((__nothrow__ , __leaf__));


extern int timer_gettime (timer_t __timerid, struct itimerspec *__value)
     __attribute__ ((__nothrow__ , __leaf__));


extern int timer_getoverrun (timer_t __timerid) __attribute__ ((__nothrow__ , __leaf__));





extern int timespec_get (struct timespec *__ts, int __base)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));
# 430 "/usr/include/time.h" 3 4

# 7 "interpolation.c" 2


# 1 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petsctime.h" 1
# 10 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petsctime.h"

# 10 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petsctime.h"
extern PetscErrorCode PetscGetCPUTime(PetscLogDouble*);


extern PetscLogDouble petsc_BaseTime;
# 106 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petsctime.h"
static inline PetscErrorCode PetscTime(PetscLogDouble *v)
{
  *v = MPI_Wtime();
  return 0;
}

static inline PetscErrorCode PetscTimeSubtract(PetscLogDouble *v)
{
  *v -= MPI_Wtime();
  return 0;
}

static inline PetscErrorCode PetscTimeAdd(PetscLogDouble *v)
{
  *v += MPI_Wtime();
  return 0;
}
# 10 "interpolation.c" 2
# 1 "/sw/eb/sw/PETSc/3.6.2-intel-2019b-Python-2.7.16/include/petscdmcomposite.h" 1






extern PetscErrorCode DMCompositeCreate(MPI_Comm,DM*);
extern PetscErrorCode DMCompositeAddDM(DM,DM);
extern PetscErrorCode DMCompositeSetCoupling(DM,PetscErrorCode (*)(DM,Mat,PetscInt*,PetscInt*,PetscInt,PetscInt,PetscInt,PetscInt));
extern PetscErrorCode DMCompositeAddVecScatter(DM,VecScatter);
extern PetscErrorCode DMCompositeScatter(DM,Vec,...);
extern PetscErrorCode DMCompositeScatterArray(DM,Vec,Vec*);
extern PetscErrorCode DMCompositeGather(DM,Vec,InsertMode,...);
extern PetscErrorCode DMCompositeGatherArray(DM,Vec,InsertMode,Vec*);
extern PetscErrorCode DMCompositeGetNumberDM(DM,PetscInt*);
extern PetscErrorCode DMCompositeGetAccess(DM,Vec,...);
extern PetscErrorCode DMCompositeRestoreAccess(DM,Vec,...);
extern PetscErrorCode DMCompositeGetAccessArray(DM,Vec,PetscInt,const PetscInt*,Vec*);
extern PetscErrorCode DMCompositeRestoreAccessArray(DM,Vec,PetscInt,const PetscInt*,Vec*);
extern PetscErrorCode DMCompositeGetLocalVectors(DM,...);
extern PetscErrorCode DMCompositeGetEntries(DM,...);
extern PetscErrorCode DMCompositeGetEntriesArray(DM,DM[]);
extern PetscErrorCode DMCompositeRestoreLocalVectors(DM,...);
extern PetscErrorCode DMCompositeGetGlobalISs(DM,IS*[]);
extern PetscErrorCode DMCompositeGetLocalISs(DM,IS**);
extern PetscErrorCode DMCompositeGetISLocalToGlobalMappings(DM,ISLocalToGlobalMapping**);
# 11 "interpolation.c" 2

PetscInt np = 0;
PetscInt ti = 0;
PetscReal L_dim = 1.0,cl = 1.0;
PetscInt block_number = 1;







PetscErrorCode ParticleVectorCreate(PetscInt numParticles, UserCtx *user) {
  PetscErrorCode ierr;
  const PetscInt particleSize = (sizeof(Particle) / sizeof(PetscReal)) + 1;

  PetscPrintf(PETSC_COMM_WORLD, "ParticleVectorCreate - cp8.1 - user - %p  \n",user);

  PetscPrintf(PETSC_COMM_WORLD, "ParticleVectorCreate - cp8.2 - particleSize - %d  \n",particleSize);

    if (!user) {
        PetscPrintf(PETSC_COMM_WORLD, "ParticleVectorCreate - cp8.3F - Error: user is NULL\n");
        return 85;
    }

      PetscPrintf(PETSC_COMM_WORLD,"ParticleVectorCreate - cp8.3S - ParticleVec before creation: %p \n", &user->ParticleVec);
      PetscPrintf(PETSC_COMM_WORLD,"ParticelVectorCrate - cp8.4 - Local ParticleVec before creation: %p \n", &user->lParticleVec);

    PetscPrintf(PETSC_COMM_WORLD,"ParticleVectorCreate - cp8.5 \n");

    PetscPrintf(PETSC_COMM_WORLD,"ParticleVectorCreate - cp8.6 - No.of Particles - %d \n",numParticles);

    if (user->ParticleVec == 
# 43 "interpolation.c" 3 4
                            ((void *)0)
# 43 "interpolation.c"
                                ) {
        PetscPrintf(PETSC_COMM_WORLD, "ParticleVectorCreate - cp8.7F - ParticleVec does not exist - Creating ParticleVec...\n");
 ierr = VecCreate(PETSC_COMM_WORLD,&user->ParticleVec); do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),45,__func__,"interpolation.c",ierr,PETSC_ERROR_REPEAT," ");} while (0);
 if (ierr) {
   PetscPrintf(PETSC_COMM_WORLD, "ParticleVectorCreate - cp8.71F - Error in DMCreateGlobalVector: %d\n", ierr);

 } else {
   PetscPrintf(PETSC_COMM_WORLD, "ParticleVectorCreate - cp8.71S - Global vector created successfully.\n");
 }

        ierr = VecSetSizes(user->ParticleVec, -1, particleSize*numParticles); do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),53,__func__,"interpolation.c",ierr,PETSC_ERROR_REPEAT," ");} while (0);
        PetscPrintf(PETSC_COMM_WORLD, "ParticleVectorCreate - cp8.72 -  vector size set \n");

        ierr = VecSetFromOptions(user->ParticleVec); do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),56,__func__,"interpolation.c",ierr,PETSC_ERROR_REPEAT," ");} while (0);
        PetscPrintf(PETSC_COMM_WORLD, "ParticleVectorCreate - cp8.73 -  vector initialized to be set from options \n");

 ierr = VecSet(user->ParticleVec,0.0); do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),59,__func__,"interpolation.c",ierr,PETSC_ERROR_REPEAT," ");} while (0);
     PetscPrintf(PETSC_COMM_WORLD, "ParticleVectorCreate - cp8.74 -  vector initialized to zeros \n");

    } else {
        PetscPrintf(PETSC_COMM_WORLD, "ParticleVectorCreate - cp8.7S - ParticleVec already exists.\n");
    }

    PetscInt N;

    VecGetSize(user->ParticleVec, &N);

    PetscPrintf(PETSC_COMM_WORLD,"ParticleVectorCreate - cp8.8 - Exit - sizeof(ParticeVec) - %d \n",N);

    return(0);
}



PetscErrorCode ParticleVectorInsert(UserCtx *user, PetscInt index, Particle *particle)
{

  PetscPrintf(PETSC_COMM_WORLD, "  ParticleVectorInsert - cp9.6%d1 - user - %p \n", index,user);

  const PetscInt particleSize = sizeof(Particle) / sizeof(PetscScalar) + 1;

  PetscPrintf(PETSC_COMM_WORLD, "  ParticleVectorInsert - cp9.6%d2 - particleSize - %d \n", index+1,particleSize);

 PetscPrintf(PETSC_COMM_WORLD, "  ParticleVectorInsert - cp9.6%d3 - particle - %p \n", index+1,particle);

  PetscInt idx[particleSize];

  PetscScalar values[particleSize];

  PetscErrorCode ierr;

  PetscPrintf(PETSC_COMM_WORLD, "  ParticleVectorInsert - cp9.6%d4 - idx - %p - values - %p \n", index+1,&idx,&values);

  values[0] = (PetscScalar)particle->cell[0];
  values[1] = (PetscScalar)particle->cell[1];
  values[2] = (PetscScalar)particle->cell[2];
  values[3] = particle->loc.x;
  values[4] = particle->loc.y;
  values[5] = particle->loc.z;
  values[6] = particle->vel.x;
  values[7] = particle->vel.y;
  values[8] = particle->vel.z;
  values[9] = particle->weights.x;
  values[10] = particle->weights.y;
  values[11] = particle->weights.z;

  PetscPrintf(PETSC_COMM_WORLD, "  ParticleVectorInsert - cp9.6%d4 - values updated - sample - values[3]- %f == loc.x - %f \n", index+1, values[3],particle->loc.x);

  for(int j = 0; j< particleSize;j++){

    idx[j] = index+j;

    ierr = VecSetValues(user->ParticleVec,1,&idx[j],&values[j],INSERT_VALUES); do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),115,__func__,"interpolation.c",ierr,PETSC_ERROR_REPEAT," ");} while (0);

  }

  PetscPrintf(PETSC_COMM_WORLD, "  ParticleVectorInsert - cp9.6%d5 - Values set into ParticleVec - %p \n", index+1,&user->ParticleVec);

  return(0);
}

PetscErrorCode ParticleVectorInitialize(UserCtx *user, PetscInt numParticles) {
  PetscRandom rand;

  PetscInt particleSize = sizeof(Particle) / sizeof(PetscReal) + 1;

    PetscPrintf(PETSC_COMM_WORLD,"ParticleVectorInitialize - cp9.1 - user - %p \n",user);
    PetscPrintf(PETSC_COMM_WORLD,"ParticleVectorInitialize - cp9.2 - No.of Particles = %d \n",numParticles);
    PetscPrintf(PETSC_COMM_WORLD,"ParticleVectorInitialize - cp9.3 - ParticleSize = %d \n",particleSize);

    srand(time(
# 133 "interpolation.c" 3 4
              ((void *)0)
# 133 "interpolation.c"
                  ));


      if (!user->ParticleVec) {
        PetscPrintf(PETSC_COMM_WORLD, "ParticleVectorInitialize - cp9.4F - Error: user->ParticleVec is NULL\n");
        return 85;
     }


    PetscRandomCreate(PETSC_COMM_WORLD, &rand);
    PetscRandomSetType(rand, "rand48");
    PetscRandomSetInterval(rand,0.0 ,1.0 );
    PetscRandomSeed(rand);

    PetscPrintf(PETSC_COMM_WORLD,"ParticleVectorInitialize - cp9.4S - Random Seed - %d \n",(PetscInt)rand);


    Particle *particleData = (Particle *)calloc(numParticles,sizeof(Particle));

PetscPrintf(PETSC_COMM_WORLD,"ParticleVectorInitialize - cp9.5 - particleData  - %p \n",particleData);

    for (PetscInt i = 0; i < numParticles; i++) {

      PetscPrintf(PETSC_COMM_WORLD,"ParticleVectorInitialize - cp9.5%dB - particleData[i]: %p - i: %d - loc.x - %f, loc.y - %f, loc.z - %f \n",i+1,particleData+i,i,particleData[i].loc.x, particleData[i].loc.y, particleData[i].loc.z);

        PetscRandomGetValue(rand, &particleData[i].loc.x);
        PetscRandomGetValue(rand, &particleData[i].loc.y);
        PetscRandomGetValue(rand, &particleData[i].loc.z);

        particleData[i].vel.x = 0.0;
        particleData[i].vel.y = 0.0;
        particleData[i].vel.z = 0.0;
        particleData[i].weights.x = 0.0;
        particleData[i].weights.y = 0.0;
        particleData[i].weights.z = 0.0;
        particleData[i].cell[0] = -1;
        particleData[i].cell[1] = -1;
        particleData[i].cell[2] = -1;


      PetscPrintf(PETSC_COMM_WORLD,"ParticleVectorInitialize - cp9.51%dA - i  %d : loc.x - %f, loc.y - %f, loc.z - %f \n",i+1,i,particleData[i].loc.x, particleData[i].loc.y, particleData[i].loc.z);
    }

   if (!user) {
    PetscPrintf(PETSC_COMM_WORLD, "ParticleVectorInitialize - cp9.6F - Error: user is NULL\n");
    return 85;
    }


    for (PetscInt i = 0; i < numParticles; ++i) {

      PetscPrintf(PETSC_COMM_WORLD," ParticleVectorInitialize - cp9.6S%d - user: %p - particleData[i]: %p - i: %d \n",i+1,user,&particleData[i],i);


      ParticleVectorInsert(user,i, &particleData[i]);
    }


   PetscPrintf(PETSC_COMM_WORLD, "ParticleVectorInitialize - cp9.7 - ParticleVector Randomized \n");


    VecAssemblyBegin(user->ParticleVec);
    VecAssemblyEnd(user->ParticleVec);

 PetscPrintf(PETSC_COMM_WORLD, "ParticleVectorInitialize - cp9.8 - ParticleVector Assembled \n");


    free(particleData);
    PetscRandomDestroy(&rand);

 PetscPrintf(PETSC_COMM_WORLD, "ParticleVectorInitialize - cp9.9 - particleData and rand freed \n");

    return(0);
}




void WriteBoundingBox(UserCtx *user, BoundingBox *bbox){

  PetscInt i,j,k,rank;
  PetscErrorCode ierr;
  PetscInt gxs,gys,gzs,gxe,gye,gze;
  DMDALocalInfo info;
  Vec Coor;
  Cmpnts ***coor,min_coords,max_coords;

  PetscPrintf(PETSC_COMM_WORLD, " WriteBoundingBox - cp10.21 - user - %p \n",user);

  DM da = user->da;
  DM fda = user->fda;

  PetscPrintf(PETSC_COMM_WORLD, " WriteBoundingBox - cp10.22 - da - %p - user->da - %p \n",da,user->da);

  PetscPrintf(PETSC_COMM_WORLD, " WriteBoundingBox - cp10.23 - fda - %p - user->fda - %p \n",fda,user->fda);

  ierr = DMGetCoordinates(user->fda, &Coor); do {if (__builtin_expect(!!(ierr),0)) return PetscError(((MPI_Comm)0x44000001),229,__func__,"interpolation.c",ierr,PETSC_ERROR_REPEAT," ");} while (0);

if (Coor == 
# 231 "interpolation.c" 3 4
           ((void *)0)
# 231 "interpolation.c"
               ) {
    PetscPrintf(PETSC_COMM_WORLD, " WriteBoundingBox - cp 10.24F - Error: Coor vector is NULL.\n");
    return;
}

  PetscPrintf(PETSC_COMM_WORLD, " WriteBoundingBox - cp10.24S - coor - %p - Coor - %p \n",&Coor,&coor);

  DMDAVecGetArrayRead(user->fda,Coor,&coor);

  PetscPrintf(PETSC_COMM_WORLD, " WriteBoundingBox - cp10.25 - fda - %p - Coor - %p - coor - %p \n",user->fda,&Coor,&coor);

  DMDAGetLocalInfo(user->da, &info);

  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;


   min_coords.x = min_coords.y = min_coords.z = 1.7976931348623157e+308/4.0;
   max_coords.x = max_coords.y = max_coords.z = -1.7976931348623157e+308/4.0;


  for (k = gzs; k < gze; k++) {
     for (j = gys; j < gye ; j++) {
         for (i = gxs; i < gxe; i++) {


                Cmpnts p0 = coor[k][j][i];
                Cmpnts p1 = coor[k][j][i + 1];
                Cmpnts p2 = coor[k][j + 1][i];
                Cmpnts p3 = coor[k][j + 1][i + 1];
                Cmpnts p4 = coor[k + 1][j][i];
                Cmpnts p5 = coor[k + 1][j][i + 1];
                Cmpnts p6 = coor[k + 1][j + 1][i];
                Cmpnts p7 = coor[k + 1][j + 1][i + 1];


                min_coords.x = (((min_coords.x)<(((((((p0.x)<(p1.x)) ? (p0.x) : (p1.x)))<((((p2.x)<(p3.x)) ? (p2.x) : (p3.x)))) ? ((((p0.x)<(p1.x)) ? (p0.x) : (p1.x))) : ((((p2.x)<(p3.x)) ? (p2.x) : (p3.x)))))) ? (min_coords.x) : (((((((p0.x)<(p1.x)) ? (p0.x) : (p1.x)))<((((p2.x)<(p3.x)) ? (p2.x) : (p3.x)))) ? ((((p0.x)<(p1.x)) ? (p0.x) : (p1.x))) : ((((p2.x)<(p3.x)) ? (p2.x) : (p3.x))))));
                min_coords.y = (((min_coords.y)<(((((((p0.y)<(p1.y)) ? (p0.y) : (p1.y)))<((((p2.y)<(p3.y)) ? (p2.y) : (p3.y)))) ? ((((p0.y)<(p1.y)) ? (p0.y) : (p1.y))) : ((((p2.y)<(p3.y)) ? (p2.y) : (p3.y)))))) ? (min_coords.y) : (((((((p0.y)<(p1.y)) ? (p0.y) : (p1.y)))<((((p2.y)<(p3.y)) ? (p2.y) : (p3.y)))) ? ((((p0.y)<(p1.y)) ? (p0.y) : (p1.y))) : ((((p2.y)<(p3.y)) ? (p2.y) : (p3.y))))));
                min_coords.z = (((min_coords.z)<(((((((p0.z)<(p1.z)) ? (p0.z) : (p1.z)))<((((p2.z)<(p3.z)) ? (p2.z) : (p3.z)))) ? ((((p0.z)<(p1.z)) ? (p0.z) : (p1.z))) : ((((p2.z)<(p3.z)) ? (p2.z) : (p3.z)))))) ? (min_coords.z) : (((((((p0.z)<(p1.z)) ? (p0.z) : (p1.z)))<((((p2.z)<(p3.z)) ? (p2.z) : (p3.z)))) ? ((((p0.z)<(p1.z)) ? (p0.z) : (p1.z))) : ((((p2.z)<(p3.z)) ? (p2.z) : (p3.z))))));

                max_coords.x = (((max_coords.x)<(((((((p0.x)<(p1.x)) ? (p1.x) : (p0.x)))<((((p2.x)<(p3.x)) ? (p3.x) : (p2.x)))) ? ((((p2.x)<(p3.x)) ? (p3.x) : (p2.x))) : ((((p0.x)<(p1.x)) ? (p1.x) : (p0.x)))))) ? (((((((p0.x)<(p1.x)) ? (p1.x) : (p0.x)))<((((p2.x)<(p3.x)) ? (p3.x) : (p2.x)))) ? ((((p2.x)<(p3.x)) ? (p3.x) : (p2.x))) : ((((p0.x)<(p1.x)) ? (p1.x) : (p0.x))))) : (max_coords.x));
                max_coords.y = (((max_coords.y)<(((((((p0.y)<(p1.y)) ? (p1.y) : (p0.y)))<((((p2.y)<(p3.y)) ? (p3.y) : (p2.y)))) ? ((((p2.y)<(p3.y)) ? (p3.y) : (p2.y))) : ((((p0.y)<(p1.y)) ? (p1.y) : (p0.y)))))) ? (((((((p0.y)<(p1.y)) ? (p1.y) : (p0.y)))<((((p2.y)<(p3.y)) ? (p3.y) : (p2.y)))) ? ((((p2.y)<(p3.y)) ? (p3.y) : (p2.y))) : ((((p0.y)<(p1.y)) ? (p1.y) : (p0.y))))) : (max_coords.y));
                max_coords.z = (((max_coords.z)<(((((((p0.z)<(p1.z)) ? (p1.z) : (p0.z)))<((((p2.z)<(p3.z)) ? (p3.z) : (p2.z)))) ? ((((p2.z)<(p3.z)) ? (p3.z) : (p2.z))) : ((((p0.z)<(p1.z)) ? (p1.z) : (p0.z)))))) ? (((((((p0.z)<(p1.z)) ? (p1.z) : (p0.z)))<((((p2.z)<(p3.z)) ? (p3.z) : (p2.z)))) ? ((((p2.z)<(p3.z)) ? (p3.z) : (p2.z))) : ((((p0.z)<(p1.z)) ? (p1.z) : (p0.z))))) : (max_coords.z));
            }
        }
    }
    DMDAVecRestoreArrayRead(fda,Coor,&coor);

    bbox->max_coords = max_coords;
    bbox->min_coords = min_coords;

}

PetscBool CPUPointIntersectCheck(BoundingBox *bbox, Particle *particle){

   Cmpnts loc,min_coords,max_coords;
   loc = particle->loc;
   PetscBool Intersects = PETSC_FALSE;
   min_coords = bbox->min_coords;
   max_coords = bbox->max_coords;

    if ((loc.x >= min_coords.x && loc.x <= max_coords.x) &&
        (loc.y >= min_coords.y && loc.y <= max_coords.y) &&
        (loc.z >= min_coords.z && loc.z <= max_coords.z)) {
        Intersects = PETSC_TRUE;
    }
   return Intersects;
}




PetscReal distance_search(Cmpnts p1, Cmpnts p2, Cmpnts p3, Cmpnts p4, Cmpnts p, PetscReal *d)
{
  PetscReal xn1, yn1, zn1;
  PetscReal xc, yc, zc;

  PetscReal dx1, dy1, dz1, dx2, dy2, dz2, r;

  dx1 = p3.x - p1.x;
  dy1 = p3.y - p1.y;
  dz1 = p3.z - p1.z;

  dx2 = p4.x - p2.x;
  dy2 = p4.y - p2.y;
  dz2 = p4.z - p2.z;

  xn1 = dy1 * dz2 - dz1 * dy2;
  yn1 = - (dx1 * dz2 - dz1 * dx2);
  zn1 = dx1 * dy2 - dy1 * dx2;

  r = sqrt(xn1 * xn1 + yn1 * yn1 + zn1 * zn1);
  xn1 /= r; yn1 /= r; zn1 /= r;

  xc = 0.25 * (p1.x + p2.x + p3.x + p4.x);
  yc = 0.25 * (p1.y + p2.y + p3.y + p4.y);
  zc = 0.25 * (p1.z + p2.z + p3.z + p4.z);

  *d = (p.x - xc) * xn1 + (p.y - yc) * yn1 + (p.z - zc) * zn1;
  if ((((*d)<0) ? -(*d) : (*d))<1.e-6) *d=0.;
  return (0);
}

PetscErrorCode CellPointInterceptCalculate(Cmpnts p, Cmpnts cell[8], PetscReal d[6])
{



  distance_search(cell[0], cell[1], cell[2], cell[3], p, &(d[4]));
  distance_search(cell[4], cell[7], cell[6], cell[5], p, &(d[5]));


  distance_search(cell[0], cell[4], cell[5], cell[1], p, &(d[2]));
  distance_search(cell[3], cell[2], cell[6], cell[7], p, &(d[3]));


  distance_search(cell[0], cell[3], cell[7], cell[4], p, &(d[0]));
  distance_search(cell[1], cell[5], cell[6], cell[2], p, &(d[1]));

  return(0);
}

PetscErrorCode CellValuesGet(Cmpnts ***coor, PetscInt idx, PetscInt idy, PetscInt idz, Cmpnts cell[8])
{
  cell[0] = coor[idz][idy][idx];
  cell[1] = coor[idz][idy][idx+1];
  cell[2] = coor[idz+1][idy][idx+1];
  cell[3] = coor[idz+1][idy][idx];
  cell[4] = coor[idz][idy+1][idx];
  cell[5] = coor[idz][idy+1][idx+1];
  cell[6] = coor[idz+1][idy+1][idx+1];
  cell[7] = coor[idz+1][idy+1][idx];

  return(0);
}

PetscBool CellIntersectCheck(PetscReal d[6])
{
  PetscBool Intersects = PETSC_FALSE;
  for(int i = 0; i<6; i++){
    if(d[i] <= 0.0) {
       Intersects = PETSC_TRUE;
       break;
    }
  }
  return Intersects;
}


PetscErrorCode InterpolationWeightsCalculate(Cmpnts a, PetscReal d[6])
{
  a.x = d[0]/(d[0]+d[1]);
  a.y = d[2]/(d[2]+d[3]);
  a.z = d[4]/(d[4]+d[5]);
  return(0);
}

PetscErrorCode CellIndexUpdate(PetscReal *d, PetscInt *idx, PetscInt *idy, PetscInt *idz)
{

  if(d[0] || d[1] <0){
    if(d[0]<0) idx--;
    else idx++;
    }


  if(d[2] || d[3] <0){
    if(d[2]<0) idy--;
    else idy++;
    }


  if(d[4] || d[5] <0){
    if(d[4]<0) idz--;
    else idz++;
    }
  return(0);
}



PetscErrorCode WalkingSearch(UserCtx *user, Particle *particle)
{
  PetscInt i,j,k;
  PetscInt idx,idy,idz;
  PetscInt xs,ys,zs,xe,ye,ze;
  PetscInt mz, my, mx;
  PetscInt lxs, lxe, lys, lye, lzs, lze;
  DM da = user->da;
  DM fda = user->fda;
  DMDALocalInfo info;
  Vec Coor;
  Cmpnts ***coor, p, cell[8];
  PetscReal d[6];
  PetscBool Cell_found = PETSC_FALSE;


  DMGetCoordinatesLocal(fda, &Coor);
  DMDAVecGetArrayRead(fda,Coor,&coor);

  DMDAGetLocalInfo(da, &info);
  mx = info.mx; my = info.my; mz = info.mz;
  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;

  lxs = xs-1; lxe = xe+1;
  lys = ys-1; lye = ye+1;
  lzs = zs-1; lze = ze+1;

  if (xs==0) lxs = xs;
  if (ys==0) lys = ys;
  if (zs==0) lzs = zs;

  if (xe==mx) lxe=xe;
  if (ye==my) lye=ye;
  if (ze==mz) lze=ze;

  idx = lxs;
  idy = lys;
  idz = lzs;
  p = particle->loc;



  while (Cell_found== PETSC_FALSE){
    CellValuesGet(coor,idx,idy,idz,cell);
    CellPointInterceptCalculate(p,cell,d);
    Cell_found = CellIntersectCheck(d);
    if(!Cell_found) CellIndexUpdate(&d,&idx,&idy,&idz);
  }

  InterpolationWeightsCalculate(particle->weights,d);

  particle->cell[0] = idx; particle->cell[1] = idy; particle->cell[2] = idz;

  DMDAVecRestoreArrayRead(fda, Coor, &coor);

  return (0);
}



PetscErrorCode ParticlesLocate(UserCtx *user, PetscInt numParticles) {

    Particle *particles;

    PetscPrintf(PETSC_COMM_WORLD, "ParticlesLocate - cp10.1 - user - %p \n",user);

    PetscInt particleSize = sizeof(Particle) / sizeof(PetscReal) + 1;

    PetscPrintf(PETSC_COMM_WORLD, "ParticlesLocate - cp10.2 - particleSize - %d - numParticles - %d \n",particleSize,numParticles);

    BoundingBox bbox;


    WriteBoundingBox(user,&bbox);


    VecGetArray(user->ParticleVec, (PetscScalar**)&particles);


    for (PetscInt i = 0; i < numParticles; ++i) {
        if (CPUPointIntersectCheck(&bbox, &particles[i])) {
            WalkingSearch(user, &particles[i]);
        }
    }


    VecRestoreArray(user->ParticleVec, (PetscScalar**)&particles);

    return(0);
}


PetscErrorCode InterpolationCoefficientsCalculate(Cmpnts a,PetscReal coeffs[8])
{
  PetscReal a1,a2,a3;
  a1 = a.x;
  a2 = a.y;
  a3 = a.z;

  coeffs[0] = a1*a2*a3;
  coeffs[1] = (a1-1.0)*a2*a3;
  coeffs[2] = a1*(a2-1.0)*a3;
  coeffs[3] = (a1-1.0)*(a2-1.0)*a3;
  coeffs[4] = a1*a2*(a3-1.0);
  coeffs[5] = (a1-1.0)*a2*(a3-1.0);
  coeffs[6] = a1*(a2-1.0)*a3;
  coeffs[7] = a1*(a2-1.0)*(a3-1.0);
  coeffs[8] = (a1-1.0)*(a2-1.0)*(a3-1.0);

  return(0);
}


PetscErrorCode InterpolationMatrixInitialize(PetscInt numParticles)
{

  PetscInt rows, cols;
  Mat InterpMat;

  PetscReal coeffs[8];

  rows = numParticles;
  cols = 8 * numParticles;


  MatCreate(PETSC_COMM_WORLD, InterpMat);
  MatSetSizes(InterpMat, -1,-1,rows,cols);
  MatSetType(InterpMat, "mpiaij");
  MatSetUp(InterpMat);
  MatZeroEntries(InterpMat);

  return(0);
}


PetscErrorCode InterpolationMatrixPopulate(Mat InterpMat, Vec particleVec,PetscInt numParticles)
{
  Particle *particles;
  PetscReal coeffs[8];
  PetscInt row;
  PetscInt colIndices[8];


  VecGetArray(particleVec, (PetscScalar**)&particles);

  for(int p = 0; p < numParticles; p++){

      row = p;

      for(int q = 0; q<8; q++){

          colIndices[q] = p*8 + q;

      }

      InterpolationCoefficientsCalculate(particles[p].weights,coeffs);

      MatSetValues(InterpMat,1, &row, 8, colIndices,coeffs,INSERT_VALUES);

  }

  VecRestoreArray(particleVec, (PetscScalar**)&particles);


  MatAssemblyBegin(InterpMat, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(InterpMat, MAT_FINAL_ASSEMBLY);

  return(0);
}


PetscErrorCode InterpolationVectorCreate(UserCtx* user, Vec particleVec, PetscInt numParticles)
{
  Particle *particles;

  PetscInt Vecsize = numParticles*8, indices[8],startid;

  PetscInt CmpntsSize = sizeof(Cmpnts) / sizeof(PetscReal);

  Vec Ucat = user->Ucat, InterpVec;

  DM da = user->da, fda = user->fda;
  DMDALocalInfo info = user->info;

  Cmpnts ***ucat,host_vel[8];

  VecCreate(((MPI_Comm)0x44000000), &InterpVec);
  VecSetSizes(InterpVec, -1,Vecsize*CmpntsSize);
  VecSetFromOptions(InterpVec);

  VecSet(InterpVec, 0.0);


  DMDAVecGetArray(fda,Ucat,&ucat);


  VecGetArray(particleVec, (PetscScalar**)&particles);

  for(int i = 0 ; i < numParticles; i++){

    startid = i;

    for( int j = 0; j < 8; j++){

        indices[j] = startid + j;

    }

      CellValuesGet(ucat,particles[i].cell[0], particles[i].cell[1], particles[i].cell[2],host_vel);
      VecSetValues(InterpVec,8,indices, host_vel,INSERT_VALUES);
  }


  DMDAVecRestoreArray(fda, Ucat, &ucat);

  VecRestoreArray(particleVec, (PetscScalar**)&particles);

  return(0);
}

PetscErrorCode Interpolate(Mat InterpMat, Vec InterpVec, PetscInt numParticles)
{

  Vec R;
  PetscInt Vecsize = numParticles;

  PetscInt CmpntsSize = sizeof(Cmpnts)/sizeof(PetscReal);


  VecCreate(((MPI_Comm)0x44000000),&R);
  VecSetSizes(R,-1,Vecsize*CmpntsSize);
  VecSetFromOptions(R);
  VecSet(R,0.0);

  MatMult(InterpMat, InterpVec, R);

  return(0);

}

PetscErrorCode ParticleVelocityUpdate(Vec ParticleVec, Vec R, PetscInt numParticles)
{

  Particle *particles;
  Cmpnts *resultant;
  PetscInt particleSize = sizeof(Particle) / sizeof(PetscReal);

  VecGetArray(ParticleVec, (PetscScalar**)&particles);
  VecGetArray(R, (PetscScalar**)&resultant);

  for(PetscInt i = 0; i<numParticles; i++){

    particles[i].vel.x = resultant[i].x;
    particles[i].vel.y = resultant[i].y;
    particles[i].vel.z = resultant[i].z;
  }


  VecRestoreArray(ParticleVec, (PetscScalar**)&particles);
  VecRestoreArray(R, (PetscScalar**)&resultant);

  return(0);
}

PetscErrorCode InterpolationMatrixDestroy(Mat InterpMat)
{

  MatDestroy(&InterpMat);

  return(0);
}


PetscErrorCode InterpolationVectorDestroy(Vec InterpVec, Vec R)
{

  VecDestroy(&InterpVec);
  VecDestroy(&R);

  return(0);

}



PetscErrorCode Ucat_Binary_Input(UserCtx *user)
{
  PetscViewer viewer;
  char filen[90];
  PetscInt bi=user->_this;

  sprintf(filen, "results/ufield%5.5d_%1.1d.dat", ti, user->_this);

  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filen, FILE_MODE_READ, &viewer);

  PetscPrintf(PETSC_COMM_WORLD, "Ucat_Binary_Input - cp7.1 - user - %p  \n",user);

  PetscInt N;

  VecGetSize(user->Ucat, &N);
  PetscPrintf(PETSC_COMM_WORLD, "Ucat_Binary_Input - cp7.2  SizeOf(ucat) - %d \n", N);
  VecLoad((user->Ucat),viewer);

  PetscViewerDestroy(&viewer);

  PetscBarrier(
# 720 "interpolation.c" 3 4
              ((void *)0)
# 720 "interpolation.c"
                        );

  return 0;

}

PetscErrorCode ReadCoordinates(UserCtx *user) {

  Cmpnts ***coor;

  Vec Coor;
  PetscInt bi, i, j, k, rank, IM, JM, KM;
  PetscReal *gc;
  FILE *fd;
  PetscReal d0 = 1.;
  PetscInt generate_grid=0, grid1d=0, nblk=block_number;

  PetscOptionsInsertFile(PETSC_COMM_WORLD, "control.dat", PETSC_TRUE);

  PetscPrintf(PETSC_COMM_WORLD, "ReadCoordinates - cp5.1 - user - %p \n", user);

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  PetscPrintf(PETSC_COMM_WORLD, "ReadCoordinates - cp5.2 - rank - %d \n", rank);

  PetscOptionsGetInt(
# 745 "interpolation.c" 3 4
                    ((void *)0)
# 745 "interpolation.c"
                              , "-grid", &generate_grid, 
# 745 "interpolation.c" 3 4
                                                         ((void *)0)
# 745 "interpolation.c"
                                                                   );

  PetscReal cl = 1.;
  PetscOptionsGetReal(
# 748 "interpolation.c" 3 4
                     ((void *)0)
# 748 "interpolation.c"
                               , "-chact_leng", &cl, 
# 748 "interpolation.c" 3 4
                                                     ((void *)0)
# 748 "interpolation.c"
                                                               );

  PetscReal L_x,L_y,L_z;


  if (generate_grid) {
    PetscOptionsGetReal(
# 754 "interpolation.c" 3 4
                       ((void *)0)
# 754 "interpolation.c"
                                 , "-L_x", &L_x, 
# 754 "interpolation.c" 3 4
                                                 ((void *)0)
# 754 "interpolation.c"
                                                           );
    PetscOptionsGetReal(
# 755 "interpolation.c" 3 4
                       ((void *)0)
# 755 "interpolation.c"
                                 , "-L_y", &L_y, 
# 755 "interpolation.c" 3 4
                                                 ((void *)0)
# 755 "interpolation.c"
                                                           );
    PetscOptionsGetReal(
# 756 "interpolation.c" 3 4
                       ((void *)0)
# 756 "interpolation.c"
                                 , "-L_z", &L_z, 
# 756 "interpolation.c" 3 4
                                                 ((void *)0)
# 756 "interpolation.c"
                                                           );

    PetscPrintf(PETSC_COMM_WORLD,"ReadCoordinates - cp5.3 - ksi eta zeta -  %le %le %le \n",L_x,L_y,L_z);

  } else {
    if (!rank) {
      fd = fopen("grid.dat", "r");
      fscanf(fd, "%i\n", &block_number);
      MPI_Bcast(&block_number, 1, ((MPI_Datatype)0x4c000405), 0, PETSC_COMM_WORLD);
    }
    else {
      MPI_Bcast(&block_number, 1, ((MPI_Datatype)0x4c000405), 0, PETSC_COMM_WORLD);
    }
  }

  PetscPrintf(PETSC_COMM_WORLD, "ReadCoordinates - cp5.4 - rank - %d \n", rank);

  PetscInt imm[block_number], kmm[block_number], jmm[block_number];
  if (generate_grid) {
    PetscOptionsGetIntArray(
# 775 "interpolation.c" 3 4
                           ((void *)0)
# 775 "interpolation.c"
                                     , "-im", imm, &nblk, 
# 775 "interpolation.c" 3 4
                                                          ((void *)0)
# 775 "interpolation.c"
                                                                    );
    PetscOptionsGetIntArray(
# 776 "interpolation.c" 3 4
                           ((void *)0)
# 776 "interpolation.c"
                                     , "-jm", jmm, &nblk, 
# 776 "interpolation.c" 3 4
                                                          ((void *)0)
# 776 "interpolation.c"
                                                                    );
    PetscOptionsGetIntArray(
# 777 "interpolation.c" 3 4
                           ((void *)0)
# 777 "interpolation.c"
                                     , "-km", kmm, &nblk, 
# 777 "interpolation.c" 3 4
                                                          ((void *)0)
# 777 "interpolation.c"
                                                                    );
  }

  PetscPrintf(PETSC_COMM_WORLD, "ReadCoordinates - cp5.5 -  imm[0] - %d; jmm[0] - %d; kmm[0] - %d  \n",imm[0],jmm[0],kmm[0]);

  for (bi=0; bi<block_number; bi++) {




    if (!rank) {
      if (!generate_grid)
 fscanf(fd, "%i %i %i\n", &(user[bi].IM), &(user[bi].JM), &(user[bi].KM));

      else {
 user[bi].IM=imm[bi];
 user[bi].JM=jmm[bi];
 user[bi].KM=kmm[bi];
      }
      IM = user[bi].IM; JM = user[bi].JM; KM = user[bi].KM;

      MPI_Bcast(&(user[bi].IM), 1, ((MPI_Datatype)0x4c000405), 0, PETSC_COMM_WORLD);
      MPI_Bcast(&(user[bi].JM), 1, ((MPI_Datatype)0x4c000405), 0, PETSC_COMM_WORLD);
      MPI_Bcast(&(user[bi].KM), 1, ((MPI_Datatype)0x4c000405), 0, PETSC_COMM_WORLD);
    }
    else {
      MPI_Bcast(&(user[bi].IM), 1, ((MPI_Datatype)0x4c000405), 0, PETSC_COMM_WORLD);
      MPI_Bcast(&(user[bi].JM), 1, ((MPI_Datatype)0x4c000405), 0, PETSC_COMM_WORLD);
      MPI_Bcast(&(user[bi].KM), 1, ((MPI_Datatype)0x4c000405), 0, PETSC_COMM_WORLD);

      IM = user[bi].IM; JM = user[bi].JM; KM = user[bi].KM;
    }

    PetscPrintf(PETSC_COMM_WORLD,"ReadCoordinates - cp5.6 - IM - %d; JM - %d; KM - %d \n",IM,JM,KM);


    PetscPrintf(PETSC_COMM_WORLD,"ReadCoordinates - cp5.61 - user[bi].da - %p;\n",&(user[bi].da));

    DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_BOX,
        user[bi].IM+1, user[bi].JM+1, user[bi].KM+1, 1,1,
   -1, 1, 2, 
# 817 "interpolation.c" 3 4
                      ((void *)0)
# 817 "interpolation.c"
                                , 
# 817 "interpolation.c" 3 4
                                  ((void *)0)
# 817 "interpolation.c"
                                            , 
# 817 "interpolation.c" 3 4
                                              ((void *)0)
# 817 "interpolation.c"
                                                        ,&(user[bi].da));

    PetscPrintf(PETSC_COMM_WORLD,"ReadCoordinates - cp5.7 - user[bi].da - %p;\n",&(user[bi].da));


    DMDASetUniformCoordinates(user[bi].da, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);

    PetscPrintf(PETSC_COMM_WORLD,"ReadCoordinates - cp5.71 - user[bi].fda - %p;\n",&(user[bi].fda));

    DMGetCoordinateDM(user[bi].da, &(user[bi].fda));

    PetscPrintf(PETSC_COMM_WORLD,"ReadCoordinates - cp5.8 - user[bi].fda - %p;\n",&(user[bi].fda));


    DMDAGetLocalInfo(user[bi].da, &(user[bi].info));

    DMDALocalInfo info = user[bi].info;
    PetscInt xs = info.xs, xe = info.xs + info.xm;
    PetscInt ys = info.ys, ye = info.ys + info.ym;
    PetscInt zs = info.zs, ze = info.zs + info.zm;
    PetscInt mx = info.mx, my = info.my, mz = info.mz;

    PetscOptionsGetInt(
# 839 "interpolation.c" 3 4
                      ((void *)0)
# 839 "interpolation.c"
                                , "-grid", &generate_grid, 
# 839 "interpolation.c" 3 4
                                                           ((void *)0)
# 839 "interpolation.c"
                                                                     );
    PetscOptionsGetInt(
# 840 "interpolation.c" 3 4
                      ((void *)0)
# 840 "interpolation.c"
                                , "-grid1d", &grid1d, 
# 840 "interpolation.c" 3 4
                                                      ((void *)0)
# 840 "interpolation.c"
                                                                );
    if (grid1d) ((*PetscTrMalloc)(((IM+JM+KM)*sizeof(PetscReal)),841,__func__,"interpolation.c",(void**)(&gc)));
    else ((*PetscTrMalloc)((3*(IM*JM*KM)*sizeof(PetscReal)),842,__func__,"interpolation.c",(void**)(&gc)));
    DMGetCoordinatesLocal(user[bi].da, &Coor);
    DMDAVecGetArray(user[bi].fda, Coor, &coor);

    if (!rank) {
      if (grid1d) {
 PetscReal xx;

 for (i=0; i<IM; i++)
   fscanf(fd, "%le %le %le\n",&gc[i],&xx,&xx);

 for (j=0; j<JM; j++)
   fscanf(fd, "%le %le %le\n",&xx,&gc[IM+j],&xx);

 for (i=0; i<KM; i++)
   fscanf(fd, "%le %le %le\n",&xx,&xx,&gc[IM+JM+i]);

 MPI_Bcast(gc, (IM+JM+KM), ((MPI_Datatype)0x4c00080b), 0, PETSC_COMM_WORLD);

 for (k=zs; k<ze; k++) {
   for (j=ys; j<ye; j++) {
     for (i=xs; i<xe; i++) {
       if (k<KM && j<JM && i<IM) {
  coor[k][j][i].x = *(gc + i)/cl*L_dim;
  coor[k][j][i].y = *(gc + IM + j)/cl*L_dim;
  coor[k][j][i].z = *(gc + IM + JM + k)/cl*L_dim;
       }
     }
   }
 }

      } else {
 for (k=0; k<KM; k++) {
   for (j=0; j<JM; j++) {
     for (i=0; i<IM; i++) {
       if (!generate_grid)
  fscanf(fd, "%le", gc + (k*(JM*IM) + j * IM + i)*3);
       else
  *(gc+(k*JM*IM + j*IM + i)*3) = L_x/(IM-1.) * i;
     }
   }
 }

 for (k=0; k<KM; k++) {
   for (j=0; j<JM; j++) {
     for (i=0; i<IM; i++) {
       if (!generate_grid)
  fscanf(fd, "%le", gc + (k*(JM*IM) + j * IM + i)*3 + 1);
       else
  *(gc+(k*JM*IM + j*IM + i)*3+1) = L_y/(JM-1.) * j;
     }
   }
 }

 for (k=0; k<KM; k++) {
   for (j=0; j<JM; j++) {
     for (i=0; i<IM; i++) {
       if (!generate_grid)
  fscanf(fd, "%le", gc + (k*(JM*IM) + j * IM + i)*3 + 2);
     else
       *(gc+(k*JM*IM + j*IM + i)*3+2) = L_z/(KM-1.) * k;
     }
   }
 }

 MPI_Bcast(gc, 3*(IM*JM*KM), ((MPI_Datatype)0x4c00080b), 0, PETSC_COMM_WORLD);

 for (k=zs; k<ze; k++) {
   for (j=ys; j<ye; j++) {
     for (i=xs; i<xe; i++) {
       if (k<KM && j<JM && i<IM) {
  coor[k][j][i].x = *(gc + (k * (IM*JM) + j * IM + i) * 3 )/cl;
  coor[k][j][i].y = *(gc + (k * (IM*JM) + j * IM + i) * 3+1)/cl;
  coor[k][j][i].z = *(gc + (k * (IM*JM) + j * IM + i) * 3+2)/cl;
       }
     }
   }
 }
      }
    }
    else {
      if (grid1d) {
 MPI_Bcast(gc, (IM+JM+KM), ((MPI_Datatype)0x4c00080b), 0, PETSC_COMM_WORLD);

 for (k=zs; k<ze; k++) {
   for (j=ys; j<ye; j++) {
     for (i=xs; i<xe; i++) {
       if (k<KM && j<JM && i<IM) {
  coor[k][j][i].x = *(gc + i)/cl*L_dim;
  coor[k][j][i].y = *(gc + IM + j)/cl*L_dim;
  coor[k][j][i].z = *(gc + IM + JM + k)/cl*L_dim;
       }
     }
   }
 }

      } else {

 MPI_Bcast(gc, 3*(IM*JM*KM), ((MPI_Datatype)0x4c00080b), 0, PETSC_COMM_WORLD);

 for (k=zs; k<ze; k++) {
   for (j=ys; j<ye; j++) {
     for (i=xs; i<xe; i++) {
       if (k<KM && j<JM && i<IM) {
  coor[k][j][i].x = *(gc + (k * (IM*JM) + j * IM + i) * 3 )/cl;
  coor[k][j][i].y = *(gc + (k * (IM*JM) + j * IM + i) * 3+1)/cl;
  coor[k][j][i].z = *(gc + (k * (IM*JM) + j * IM + i) * 3+2)/cl;
       }
     }
   }
 }
      }
    }
    ((*PetscTrFree)((void*)(gc),955,__func__,"interpolation.c") || ((gc) = 0,0));
    DMDAVecRestoreArray(user[bi].fda, Coor, &coor);

    Vec gCoor;
    DMGetCoordinates(user[bi].da, &gCoor);
    DMLocalToGlobalBegin(user[bi].fda, Coor, INSERT_VALUES, gCoor);
    DMLocalToGlobalEnd(user[bi].fda, Coor, INSERT_VALUES, gCoor);

    DMGlobalToLocalBegin(user[bi].fda, gCoor, INSERT_VALUES, Coor);
    DMGlobalToLocalEnd(user[bi].fda, gCoor, INSERT_VALUES, Coor);
  }

  if (!rank) {
    if(!generate_grid)
      fclose(fd);
  }

  for (bi=0; bi<block_number; bi++) {
    user[bi]._this = bi;
    PetscPrintf(PETSC_COMM_WORLD,"ReadCoordinates - cp5.9%d - user[bi] - %p \n",bi+1,&user[bi]);
  }

  return(0);
}
# 987 "interpolation.c"
int main(int argc, char **argv){

DM da,fda,pda;
UserCtx *user;

PetscInitialize(&argc, &argv, (char *)0, help);
PetscOptionsInsertFile(PETSC_COMM_WORLD, "control.dat", PETSC_TRUE);

PetscInt rank, bi, ibi,i;

PetscPrintf(PETSC_COMM_WORLD, "main - cp1 \n");



((*PetscTrMalloc)((block_number*sizeof(UserCtx)),1001,__func__,"interpolation.c",(void**)(&user)));

PetscPrintf(PETSC_COMM_WORLD, "main - cp2 - user - %p  \n",user);

MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

PetscPrintf(PETSC_COMM_WORLD, "main - cp3 - rank - %d  \n", rank);





PetscOptionsGetInt(
# 1013 "interpolation.c" 3 4
                  ((void *)0)
# 1013 "interpolation.c"
                            , "-ti", &ti, 
# 1013 "interpolation.c" 3 4
                                          ((void *)0)
# 1013 "interpolation.c"
                                                    );

PetscPrintf(PETSC_COMM_WORLD, "main - cp4 - ti - %d \n", ti);

PetscOptionsGetInt(
# 1017 "interpolation.c" 3 4
                  ((void *)0)
# 1017 "interpolation.c"
                            , "-numParticles", &np, 
# 1017 "interpolation.c" 3 4
                                                    ((void *)0)
# 1017 "interpolation.c"
                                                              );

PetscPrintf(PETSC_COMM_WORLD, "main - cp5 - np - %d  \n",np);



ReadCoordinates(user);

PetscPrintf(PETSC_COMM_WORLD, "main - cp6 \n");



for(bi = 0; bi < block_number; bi ++ ){

  DMCreateGlobalVector(user[bi].fda, &user[bi].Ucat);

}

PetscPrintf(PETSC_COMM_WORLD, "main - cp7 \n");



Ucat_Binary_Input(user);

PetscPrintf(PETSC_COMM_WORLD, "main - cp8 \n");



PetscPrintf(PETSC_COMM_WORLD, "main - cp 8A - np - %d  \n",np);

PetscPrintf(PETSC_COMM_WORLD, "main - cp 8B - user - %p  \n",user);

ParticleVectorCreate(np, user);

PetscPrintf(PETSC_COMM_WORLD, "main - cp9 \n");



ParticleVectorInitialize(user,np);

PetscPrintf(PETSC_COMM_WORLD, "main - cp10 \n");



ParticlesLocate(user,np);

PetscPrintf(PETSC_COMM_WORLD, "main - cp10 \n");



((*PetscTrFree)((void*)(user),1067,__func__,"interpolation.c") || ((user) = 0,0));

PetscFinalize();

return 0;
}
