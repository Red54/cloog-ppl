
   /**-------------------------------------------------------------------**
    **                               CLooG                               **
    **-------------------------------------------------------------------**
    **                             domain.h                              **
    **-------------------------------------------------------------------**
    **                  First version: october 28th 2001                 **
    **-------------------------------------------------------------------**/


/******************************************************************************
 *               CLooG : the Chunky Loop Generator (experimental)             *
 ******************************************************************************
 *                                                                            *
 * Copyright (C) 2001-2005 Cedric Bastoul                                     *
 *                                                                            *
 * This is free software; you can redistribute it and/or modify it under the  *
 * terms of the GNU General Public License as published by the Free Software  *
 * Foundation; either version 2 of the License, or (at your option) any later *
 * version.                                                                   *
 *                                                                            *
 * This software is distributed in the hope that it will be useful, but       *
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY *
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License   *
 * for more details.                                                          *
 *                                                                            *
 * You should have received a copy of the GNU General Public License along    *
 * with software; if not, write to the Free Software Foundation, Inc.,        *
 * 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA                     *
 *                                                                            *
 * CLooG, the Chunky Loop Generator                                           *
 * Written by Cedric Bastoul, Cedric.Bastoul@inria.fr                         *
 *                                                                            *
 ******************************************************************************/


#ifndef CLOOG_DOMAIN_H
#define CLOOG_DOMAIN_H
#if defined(__cplusplus)
extern "C" 
  {
#endif 


/* The Polyhedron structure comes directly from PolyLib (defined in
 * polylib/types.h) here is how it looks like (at least in PolyLib 5.20.0
 * version).
 *
 * typedef struct polyhedron { 
 *   unsigned Dimension,      // Dimension number (NbColumns-2 in Matrix).
 *            NbConstraints,  // Number of constraints (NbRows in Matrix).
 *            NbRays,         // Number of rays in dual representation.
 *            NbEq,           // Number of vertices (?).
 *            NbBid;          // Number of extremal rays (?).
 *   Value **Constraint;      // The pointers to rows in matrix representation.
 *   Value **Ray;             // The pointers to rays in dual representation.
 *   Value *p_Init;           // The whole data, consecutive in memory.
 *   int p_Init_size;         // To clear values in GMP mode.
 *   struct polyhedron *next; // Pointer to next component of the union.
 * } Polyhedron;
 */ 


/**
 * CloogDomain structure:
 * this structure contains a polyhedron in the PolyLib shape and the number of
 * active references to this structure. Because CLooG uses many copies of
 * domains there is no need to actually copy these domains but just to return
 * a pointer to them and to increment the number of active references. Each time
 * a CloogDomain will be freed, we will decrement the active reference counter
 * and actually free it if its value is zero.
 */
struct cloogdomain
{ Polyhedron * _polyhedron ;      /**< The polyhedral domain. */
  int _references ;               /**< Number of references to this structure. */
} ;
typedef struct cloogdomain CloogDomain ;


/**
 * CloogDomainList structure:
 * this structure reprensents a node of a linked list of CloogDomain structures.
 */
struct cloogdomainlist
{ CloogDomain * _domain ;         /**< An element of the list. */
  struct cloogdomainlist * _next ;/**< Pointer to the next element of the list.*/
} ;
typedef struct cloogdomainlist CloogDomainList ;

static inline CloogDomain *cloog_domain (CloogDomainList *l)
{
  return l->_domain;
}

static inline void cloog_set_domain (CloogDomainList *l, CloogDomain *d)
{
  l->_domain = d;
}

static inline CloogDomainList *cloog_next_domain (CloogDomainList *l)
{
  return l->_next;
}

static inline void cloog_set_next_domain (CloogDomainList *l, CloogDomainList *n)
{
  l->_next = n;
}


/******************************************************************************
 *                         Polyhedral Library interface                       *
 ******************************************************************************/
void          cloog_domain_print(FILE *, CloogDomain *) ;
void          cloog_domain_free(CloogDomain *) ;
CloogDomain * cloog_domain_copy(CloogDomain *) ;
CloogDomain * cloog_domain_convex(CloogDomain *) ;
CloogDomain * cloog_domain_simple_convex(CloogDomain * domain, int nb_par);
CloogDomain * cloog_domain_simplify(CloogDomain *, CloogDomain *) ;
CloogDomain * cloog_domain_union(CloogDomain *, CloogDomain *) ;
CloogDomain * cloog_domain_intersection(CloogDomain *, CloogDomain *) ;
CloogDomain * cloog_domain_difference(CloogDomain *, CloogDomain *) ;
CloogDomain * cloog_domain_addconstraints(CloogDomain *, CloogDomain *) ;
void          cloog_domain_sort(CloogDomain **, unsigned, unsigned, unsigned, int *);
CloogDomain * cloog_domain_empty(int) ;
int cloog_domain_isconvex (CloogDomain *);
unsigned cloog_domain_dim (CloogDomain *);
unsigned cloog_domain_nb_polyhedra (CloogDomain *);
void cloog_domain_print_polyhedra (FILE *, CloogDomain *);


/******************************************************************************
 *                          Structure display function                        *
 ******************************************************************************/
void cloog_domain_print_structure(FILE *, CloogDomain *, int) ;
void cloog_domain_list_print(FILE *, CloogDomainList *) ;


/******************************************************************************
 *                         Memory deallocation function                       *
 ******************************************************************************/
void cloog_domain_list_free(CloogDomainList *) ;


/*+****************************************************************************
 *                               Reading function                             *
 ******************************************************************************/
CloogDomain * cloog_domain_read(FILE *) ;
CloogDomain * cloog_domain_union_read(FILE *) ;
CloogDomainList * cloog_domain_list_read(FILE *) ;


/******************************************************************************
 *                            Processing functions                            *
 ******************************************************************************/
int           cloog_domain_isempty(CloogDomain *) ;
CloogDomain * cloog_domain_project(CloogDomain *, int, int) ;
CloogDomain * cloog_domain_extend(CloogDomain *, int, int) ;
int           cloog_domain_never_integral(CloogDomain *) ;
void          cloog_domain_stride(CloogDomain *, int, int, Value *, Value *) ;
int           cloog_domain_integral_lowerbound(CloogDomain *, int, Value *) ;
void          cloog_domain_lowerbound_update(CloogDomain *, int, Value) ;
int           cloog_domain_lazy_disjoint(CloogDomain *, CloogDomain *) ;
int           cloog_domain_lazy_equal(CloogDomain *, CloogDomain *) ;
int           cloog_domain_lazy_block(CloogDomain *, CloogDomain *,
                                      CloogDomainList *, int) ;
int           cloog_domain_lazy_isscalar(CloogDomain *, int) ;
int           cloog_domain_list_lazy_same(CloogDomainList *) ;
void          cloog_domain_scalar(CloogDomain *, int, Value *) ;
CloogDomain * cloog_domain_cut_first(CloogDomain *) ;
CloogDomain * cloog_domain_erase_dimension(CloogDomain *, int) ;

#if defined(__cplusplus)
  }
#endif 
#endif /* define _H */
