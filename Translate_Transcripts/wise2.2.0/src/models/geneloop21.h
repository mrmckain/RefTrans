#ifndef DYNAMITEgeneloop21HEADERFILE
#define DYNAMITEgeneloop21HEADERFILE
#ifdef _cplusplus
extern "C" {
#endif

#include "dyna.h"
#include "genewisemodel.h"



struct Wise2_GeneLoop21 {  
    int dynamite_hard_link;  
    BaseMatrix * basematrix;     
    int leni;    
    int lenj;    
    GeneWiseScore * query;   
    ComplexSequence* target;     
    GeneParser21Score * gp21;    
    RandomCodonScore * rndcodon;     
    RandomModelDNAScore * rndbase;   
    } ;  
/* GeneLoop21 defined */ 
#ifndef DYNAMITE_DEFINED_GeneLoop21
typedef struct Wise2_GeneLoop21 Wise2_GeneLoop21;
#define GeneLoop21 Wise2_GeneLoop21
#define DYNAMITE_DEFINED_GeneLoop21
#endif


#ifdef PTHREAD
struct thread_pool_holder_GeneLoop21 {  
    GeneWiseScore * query;  /* Static query data: never free'd */ 
    ComplexSequence* target;/* Target object placeholder */ 
    GenomicDB* targetdb;/* Target database object */ 
    boolean target_init; 
    GeneParser21Score * gp21;    
    RandomCodonScore * rndcodon;     
    RandomModelDNAScore * rndbase;   
    pthread_mutex_t input_lock;  
    pthread_mutex_t output_lock;     
    Hscore * out;    
    pthread_t * pool;    
    int number_of_threads;   
    boolean search_has_ended;    
    DBSearchImpl * dbsi;     
    } ;  
#endif /* PTHREAD */


    /***************************************************/
    /* Callable functions                              */
    /* These are the functions you are expected to use */
    /***************************************************/



/* Function:  search_GeneLoop21(dbsi,out,query,targetdb,gp21,rndcodon,rndbase)
 *
 * Descrip:    This function makes a database search of GeneLoop21
 *             It uses the dbsi structure to choose which implementation
 *             to use of the database searching. This way at run time you
 *             can switch between single threaded/multi-threaded or hardware
 *
 *
 * Arg:            dbsi [UNKN ] Undocumented argument [DBSearchImpl *]
 * Arg:             out [UNKN ] Undocumented argument [Hscore *]
 * Arg:           query [UNKN ] Undocumented argument [GeneWiseScore *]
 * Arg:        targetdb [UNKN ] Undocumented argument [GenomicDB*]
 * Arg:            gp21 [UNKN ] Undocumented argument [GeneParser21Score *]
 * Arg:        rndcodon [UNKN ] Undocumented argument [RandomCodonScore *]
 * Arg:         rndbase [UNKN ] Undocumented argument [RandomModelDNAScore *]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_search_GeneLoop21(DBSearchImpl * dbsi,Hscore * out,GeneWiseScore * query,GenomicDB* targetdb ,GeneParser21Score * gp21,RandomCodonScore * rndcodon,RandomModelDNAScore * rndbase);
#define search_GeneLoop21 Wise2_search_GeneLoop21


/* Function:  serial_search_GeneLoop21(out,query,targetdb,gp21,rndcodon,rndbase)
 *
 * Descrip:    This function makes a database search of GeneLoop21
 *             It is a single processor implementation
 *
 *
 * Arg:             out [UNKN ] Undocumented argument [Hscore *]
 * Arg:           query [UNKN ] Undocumented argument [GeneWiseScore *]
 * Arg:        targetdb [UNKN ] Undocumented argument [GenomicDB*]
 * Arg:            gp21 [UNKN ] Undocumented argument [GeneParser21Score *]
 * Arg:        rndcodon [UNKN ] Undocumented argument [RandomCodonScore *]
 * Arg:         rndbase [UNKN ] Undocumented argument [RandomModelDNAScore *]
 *
 * Return [UNKN ]  Undocumented return value [Search_Return_Type]
 *
 */
Search_Return_Type Wise2_serial_search_GeneLoop21(Hscore * out,GeneWiseScore * query,GenomicDB* targetdb ,GeneParser21Score * gp21,RandomCodonScore * rndcodon,RandomModelDNAScore * rndbase);
#define serial_search_GeneLoop21 Wise2_serial_search_GeneLoop21


/* Function:  PackAln_bestmemory_GeneLoop21(query,target,gp21,rndcodon,rndbase,dpenv,dpri)
 *
 * Descrip:    This function chooses the best memory set-up for the alignment
 *             using calls to basematrix, and then implements either a large
 *             or small memory model.
 *
 *             It is the best function to use if you just want an alignment
 *
 *             If you want a label alignment, you will need
 *             /convert_PackAln_to_AlnBlock_GeneLoop21
 *
 *
 * Arg:           query [UNKN ] query data structure [GeneWiseScore *]
 * Arg:          target [UNKN ] target data structure [ComplexSequence*]
 * Arg:            gp21 [UNKN ] Resource [GeneParser21Score *]
 * Arg:        rndcodon [UNKN ] Resource [RandomCodonScore *]
 * Arg:         rndbase [UNKN ] Resource [RandomModelDNAScore *]
 * Arg:           dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 * Arg:            dpri [UNKN ] Undocumented argument [DPRunImpl *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_bestmemory_GeneLoop21(GeneWiseScore * query,ComplexSequence* target ,GeneParser21Score * gp21,RandomCodonScore * rndcodon,RandomModelDNAScore * rndbase,DPEnvelope * dpenv,DPRunImpl * dpri);
#define PackAln_bestmemory_GeneLoop21 Wise2_PackAln_bestmemory_GeneLoop21


/* Function:  allocate_Expl_GeneLoop21(query,target,gp21,rndcodon,rndbase)
 *
 * Descrip:    This function allocates the GeneLoop21 structure
 *             and the basematrix area for explicit memory implementations
 *             It calls /allocate_GeneLoop21_only
 *
 *
 * Arg:           query [UNKN ] query data structure [GeneWiseScore *]
 * Arg:          target [UNKN ] target data structure [ComplexSequence*]
 * Arg:            gp21 [UNKN ] Resource [GeneParser21Score *]
 * Arg:        rndcodon [UNKN ] Resource [RandomCodonScore *]
 * Arg:         rndbase [UNKN ] Resource [RandomModelDNAScore *]
 *
 * Return [UNKN ]  Undocumented return value [GeneLoop21 *]
 *
 */
GeneLoop21 * Wise2_allocate_Expl_GeneLoop21(GeneWiseScore * query,ComplexSequence* target ,GeneParser21Score * gp21,RandomCodonScore * rndcodon,RandomModelDNAScore * rndbase);
#define allocate_Expl_GeneLoop21 Wise2_allocate_Expl_GeneLoop21


/* Function:  recalculate_PackAln_GeneLoop21(pal,mat)
 *
 * Descrip:    This function recalculates the PackAln structure produced by GeneLoop21
 *             For example, in linear space methods this is used to score them
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 * Arg:        mat [UNKN ] Undocumented argument [GeneLoop21 *]
 *
 */
void Wise2_recalculate_PackAln_GeneLoop21(PackAln * pal,GeneLoop21 * mat);
#define recalculate_PackAln_GeneLoop21 Wise2_recalculate_PackAln_GeneLoop21


/* Function:  allocate_Small_GeneLoop21(query,target,gp21,rndcodon,rndbase)
 *
 * Descrip:    This function allocates the GeneLoop21 structure
 *             and the basematrix area for a small memory implementations
 *             It calls /allocate_GeneLoop21_only
 *
 *
 * Arg:           query [UNKN ] query data structure [GeneWiseScore *]
 * Arg:          target [UNKN ] target data structure [ComplexSequence*]
 * Arg:            gp21 [UNKN ] Resource [GeneParser21Score *]
 * Arg:        rndcodon [UNKN ] Resource [RandomCodonScore *]
 * Arg:         rndbase [UNKN ] Resource [RandomModelDNAScore *]
 *
 * Return [UNKN ]  Undocumented return value [GeneLoop21 *]
 *
 */
GeneLoop21 * Wise2_allocate_Small_GeneLoop21(GeneWiseScore * query,ComplexSequence* target ,GeneParser21Score * gp21,RandomCodonScore * rndcodon,RandomModelDNAScore * rndbase);
#define allocate_Small_GeneLoop21 Wise2_allocate_Small_GeneLoop21


/* Function:  PackAln_calculate_Small_GeneLoop21(mat,dpenv)
 *
 * Descrip:    This function calculates an alignment for GeneLoop21 structure in linear space
 *             If you want only the start/end points
 *             use /AlnRangeSet_calculate_Small_GeneLoop21 
 *
 *             The function basically
 *               finds start/end points 
 *               foreach start/end point 
 *                 calls /full_dc_GeneLoop21 
 *
 *
 * Arg:          mat [UNKN ] Undocumented argument [GeneLoop21 *]
 * Arg:        dpenv [UNKN ] Undocumented argument [DPEnvelope *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_calculate_Small_GeneLoop21(GeneLoop21 * mat,DPEnvelope * dpenv);
#define PackAln_calculate_Small_GeneLoop21 Wise2_PackAln_calculate_Small_GeneLoop21


/* Function:  AlnRangeSet_calculate_Small_GeneLoop21(mat)
 *
 * Descrip:    This function calculates an alignment for GeneLoop21 structure in linear space
 *             If you want the full alignment, use /PackAln_calculate_Small_GeneLoop21 
 *             If you have already got the full alignment, but want the range set, use /AlnRangeSet_from_PackAln_GeneLoop21
 *             If you have got the small matrix but not the alignment, use /AlnRangeSet_from_GeneLoop21 
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [GeneLoop21 *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_calculate_Small_GeneLoop21(GeneLoop21 * mat);
#define AlnRangeSet_calculate_Small_GeneLoop21 Wise2_AlnRangeSet_calculate_Small_GeneLoop21


/* Function:  AlnRangeSet_from_GeneLoop21(mat)
 *
 * Descrip:    This function reads off a start/end structure
 *             for GeneLoop21 structure in linear space
 *             If you want the full alignment use
 *             /PackAln_calculate_Small_GeneLoop21 
 *             If you have not calculated the matrix use
 *             /AlnRange_calculate_Small_GeneLoop21
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [GeneLoop21 *]
 *
 * Return [UNKN ]  Undocumented return value [AlnRangeSet *]
 *
 */
AlnRangeSet * Wise2_AlnRangeSet_from_GeneLoop21(GeneLoop21 * mat);
#define AlnRangeSet_from_GeneLoop21 Wise2_AlnRangeSet_from_GeneLoop21


/* Function:  convert_PackAln_to_AlnBlock_GeneLoop21(pal)
 *
 * Descrip:    Converts a path alignment to a label alignment
 *             The label alignment is probably much more useful than the path
 *
 *
 * Arg:        pal [UNKN ] Undocumented argument [PackAln *]
 *
 * Return [UNKN ]  Undocumented return value [AlnBlock *]
 *
 */
AlnBlock * Wise2_convert_PackAln_to_AlnBlock_GeneLoop21(PackAln * pal);
#define convert_PackAln_to_AlnBlock_GeneLoop21 Wise2_convert_PackAln_to_AlnBlock_GeneLoop21


/* Function:  PackAln_read_Expl_GeneLoop21(mat)
 *
 * Descrip:    Reads off PackAln from explicit matrix structure
 *
 *
 * Arg:        mat [UNKN ] Undocumented argument [GeneLoop21 *]
 *
 * Return [UNKN ]  Undocumented return value [PackAln *]
 *
 */
PackAln * Wise2_PackAln_read_Expl_GeneLoop21(GeneLoop21 * mat);
#define PackAln_read_Expl_GeneLoop21 Wise2_PackAln_read_Expl_GeneLoop21


/* Function:  calculate_GeneLoop21(mat)
 *
 * Descrip:    This function calculates the GeneLoop21 matrix when in explicit mode
 *             To allocate the matrix use /allocate_Expl_GeneLoop21
 *
 *
 * Arg:        mat [UNKN ] GeneLoop21 which contains explicit basematrix memory [GeneLoop21 *]
 *
 * Return [UNKN ]  Undocumented return value [boolean]
 *
 */
boolean Wise2_calculate_GeneLoop21(GeneLoop21 * mat);
#define calculate_GeneLoop21 Wise2_calculate_GeneLoop21


/* Function:  GeneLoop21_alloc(void)
 *
 * Descrip:    Allocates structure: assigns defaults if given 
 *
 *
 *
 * Return [UNKN ]  Undocumented return value [GeneLoop21 *]
 *
 */
GeneLoop21 * Wise2_GeneLoop21_alloc(void);
#define GeneLoop21_alloc Wise2_GeneLoop21_alloc


/* Function:  free_GeneLoop21(obj)
 *
 * Descrip:    Free Function: removes the memory held by obj
 *             Will chain up to owned members and clear all lists
 *
 *
 * Arg:        obj [UNKN ] Object that is free'd [GeneLoop21 *]
 *
 * Return [UNKN ]  Undocumented return value [GeneLoop21 *]
 *
 */
GeneLoop21 * Wise2_free_GeneLoop21(GeneLoop21 * obj);
#define free_GeneLoop21 Wise2_free_GeneLoop21


  /* Unplaced functions */
  /* There has been no indication of the use of these functions */


    /***************************************************/
    /* Internal functions                              */
    /* you are not expected to have to call these      */
    /***************************************************/
int Wise2_score_only_GeneLoop21(GeneWiseScore * query,ComplexSequence* target ,GeneParser21Score * gp21,RandomCodonScore * rndcodon,RandomModelDNAScore * rndbase);
#define score_only_GeneLoop21 Wise2_score_only_GeneLoop21
GeneLoop21 * Wise2_allocate_GeneLoop21_only(GeneWiseScore * query,ComplexSequence* target ,GeneParser21Score * gp21,RandomCodonScore * rndcodon,RandomModelDNAScore * rndbase);
#define allocate_GeneLoop21_only Wise2_allocate_GeneLoop21_only
void Wise2_init_GeneLoop21(GeneLoop21 * mat);
#define init_GeneLoop21 Wise2_init_GeneLoop21
AlnRange * Wise2_AlnRange_build_GeneLoop21(GeneLoop21 * mat,int stopj,int stopspecstate,int * startj,int * startspecstate);
#define AlnRange_build_GeneLoop21 Wise2_AlnRange_build_GeneLoop21
boolean Wise2_read_hidden_GeneLoop21(GeneLoop21 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out);
#define read_hidden_GeneLoop21 Wise2_read_hidden_GeneLoop21
int Wise2_max_hidden_GeneLoop21(GeneLoop21 * mat,int hiddenj,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_hidden_GeneLoop21 Wise2_max_hidden_GeneLoop21
boolean Wise2_read_special_strip_GeneLoop21(GeneLoop21 * mat,int stopi,int stopj,int stopstate,int * startj,int * startstate,PackAln * out);
#define read_special_strip_GeneLoop21 Wise2_read_special_strip_GeneLoop21
int Wise2_max_special_strip_GeneLoop21(GeneLoop21 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_special_strip_GeneLoop21 Wise2_max_special_strip_GeneLoop21
int Wise2_max_matrix_to_special_GeneLoop21(GeneLoop21 * mat,int i,int j,int state,int cscore,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_matrix_to_special_GeneLoop21 Wise2_max_matrix_to_special_GeneLoop21
void Wise2_calculate_hidden_GeneLoop21(GeneLoop21 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv);
#define calculate_hidden_GeneLoop21 Wise2_calculate_hidden_GeneLoop21
void Wise2_init_hidden_GeneLoop21(GeneLoop21 * mat,int starti,int startj,int stopi,int stopj);
#define init_hidden_GeneLoop21 Wise2_init_hidden_GeneLoop21
boolean Wise2_full_dc_GeneLoop21(GeneLoop21 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,PackAln * out,int * donej,int totalj,DPEnvelope * dpenv);
#define full_dc_GeneLoop21 Wise2_full_dc_GeneLoop21
boolean Wise2_do_dc_single_pass_GeneLoop21(GeneLoop21 * mat,int starti,int startj,int startstate,int stopi,int stopj,int stopstate,DPEnvelope * dpenv,int perc_done);
#define do_dc_single_pass_GeneLoop21 Wise2_do_dc_single_pass_GeneLoop21
void Wise2_push_dc_at_merge_GeneLoop21(GeneLoop21 * mat,int starti,int stopi,int startj,int * stopj,DPEnvelope * dpenv);
#define push_dc_at_merge_GeneLoop21 Wise2_push_dc_at_merge_GeneLoop21
void Wise2_follow_on_dc_GeneLoop21(GeneLoop21 * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define follow_on_dc_GeneLoop21 Wise2_follow_on_dc_GeneLoop21
void Wise2_run_up_dc_GeneLoop21(GeneLoop21 * mat,int starti,int stopi,int startj,int stopj,DPEnvelope * dpenv,int perc_done);
#define run_up_dc_GeneLoop21 Wise2_run_up_dc_GeneLoop21
void Wise2_init_dc_GeneLoop21(GeneLoop21 * mat);
#define init_dc_GeneLoop21 Wise2_init_dc_GeneLoop21
boolean Wise2_dc_start_end_calculate_GeneLoop21(GeneLoop21 * mat,DPEnvelope * dpenv);
#define dc_start_end_calculate_GeneLoop21 Wise2_dc_start_end_calculate_GeneLoop21
int Wise2_start_end_find_end_GeneLoop21(GeneLoop21 * mat,int * endj);
#define start_end_find_end_GeneLoop21 Wise2_start_end_find_end_GeneLoop21
void Wise2_init_start_end_linear_GeneLoop21(GeneLoop21 * mat);
#define init_start_end_linear_GeneLoop21 Wise2_init_start_end_linear_GeneLoop21
AlnConvertSet * Wise2_AlnConvertSet_GeneLoop21(void);
#define AlnConvertSet_GeneLoop21 Wise2_AlnConvertSet_GeneLoop21
int Wise2_find_end_GeneLoop21(GeneLoop21 * mat,int * ri,int * rj,int * state,boolean * isspecial);
#define find_end_GeneLoop21 Wise2_find_end_GeneLoop21
void Wise2_GeneLoop21_debug_show_matrix(GeneLoop21 * mat,int starti,int stopi,int startj,int stopj,FILE * ofp);
#define GeneLoop21_debug_show_matrix Wise2_GeneLoop21_debug_show_matrix
int Wise2_max_calc_GeneLoop21(GeneLoop21 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_calc_GeneLoop21 Wise2_max_calc_GeneLoop21
int Wise2_max_calc_special_GeneLoop21(GeneLoop21 * mat,int i,int j,int state,boolean isspecial,int * reti,int * retj,int * retstate,boolean * retspecial,int * cellscore);
#define max_calc_special_GeneLoop21 Wise2_max_calc_special_GeneLoop21

#ifdef _cplusplus
}
#endif

#endif
