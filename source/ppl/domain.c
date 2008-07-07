
   /**-------------------------------------------------------------------**
    **                               CLooG                               **
    **-------------------------------------------------------------------**
    **                             domain.c                              **
    **-------------------------------------------------------------------**
    **                   First version: october 28th 2001                **
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
/* CAUTION: the english used for comments is probably the worst you ever read,
 *          please feel free to correct and improve it !
 */

# include <stdlib.h>
# include <stdio.h>
# include <ctype.h>
# include "../../include/cloog/cloog.h"
#include "matrix.h"

static int cloog_check_polyhedral_ops = 1;
static int cloog_return_ppl_result = 0;
static int cloog_print_debug = 0;

static CloogDomain *
print_result (char *s, CloogDomain *d)
{
  if (cloog_print_debug)
    {
      fprintf (stderr, "%s \n", s);
      debug_cloog_domain (d);
    }
  return d;
}

/* Variables names for pretty printing.  */
static char wild_name[200][40];

static inline const char*
variable_output_function (ppl_dimension_type var)
{
  if (var < 40)
    return wild_name[var + 1];
  else
    return 0;
}

static inline void
error_handler (enum ppl_enum_error_code code, const char* description)
{
  fprintf (stderr, "PPL error code %d\n%s", code, description);
  exit (1);
}

void
cloog_initialize (void)
{
  sprintf (wild_name[0], "1");
  sprintf (wild_name[1], "a");
  sprintf (wild_name[2], "b");
  sprintf (wild_name[3], "c");
  sprintf (wild_name[4], "d");
  sprintf (wild_name[5], "e");
  sprintf (wild_name[6], "f");
  sprintf (wild_name[7], "g");
  sprintf (wild_name[8], "h");
  sprintf (wild_name[9], "i");
  sprintf (wild_name[10], "j");
  sprintf (wild_name[11], "k");
  sprintf (wild_name[12], "l");
  sprintf (wild_name[13], "m");
  sprintf (wild_name[14], "n");
  sprintf (wild_name[15], "o");
  sprintf (wild_name[16], "p");
  sprintf (wild_name[17], "q");
  sprintf (wild_name[18], "r");
  sprintf (wild_name[19], "s");
  sprintf (wild_name[20], "t");
  sprintf (wild_name[21], "alpha");
  sprintf (wild_name[22], "beta");
  sprintf (wild_name[23], "gamma");
  sprintf (wild_name[24], "delta");
  sprintf (wild_name[25], "tau");
  sprintf (wild_name[26], "sigma");
  sprintf (wild_name[27], "chi");
  sprintf (wild_name[28], "omega");
  sprintf (wild_name[29], "pi");
  sprintf (wild_name[30], "ni");
  sprintf (wild_name[31], "Alpha");
  sprintf (wild_name[32], "Beta");
  sprintf (wild_name[33], "Gamma");
  sprintf (wild_name[34], "Delta");
  sprintf (wild_name[35], "Tau");
  sprintf (wild_name[36], "Sigma");
  sprintf (wild_name[37], "Chi");
  sprintf (wild_name[38], "Omega");
  sprintf (wild_name[39], "xxx");

  if (ppl_initialize() < 0)
    {
      fprintf (stderr, "Cannot initialize the Parma Polyhedra Library.\n");
      exit (1);
    }

  if (ppl_set_error_handler (error_handler) < 0)
    {
      fprintf (stderr, "Cannot install the custom error handler.\n");
      exit (1);
    }

  if (ppl_io_set_variable_output_function (variable_output_function) < 0)
    {
      fprintf (stderr, "Cannot install the PPL custom variable output function. \n");
      exit (1);
    }
}

static inline Polyhedron *
u2p (ppl_polyhedra_union * upol)
{
  Polyhedron *res = Polyhedron_Copy (cloog_upol_polyhedron (upol));
  Polyhedron *p = res;

  while (upol)
    {
      ppl_polyhedra_union *next = cloog_upol_next (upol);
      Polyhedron *n;

      if (next)
	n = Polyhedron_Copy (cloog_upol_polyhedron (next));
      else
	n = NULL;

      p->next = n;
      p = n;
      upol = next;
    }

  return res;
}

static inline Polyhedron *
d2p (CloogDomain * d)
{
  return u2p (cloog_domain_upol (d));
}


static inline ppl_polyhedra_union *
p2u (Polyhedron * p)
{
  ppl_polyhedra_union *u = cloog_new_upol (p);
  ppl_polyhedra_union *res = u;

  while (p)
    {
      Polyhedron *next = p->next;
      ppl_polyhedra_union *n;

      if (next)
	n = cloog_new_upol (next);
      else
	n = NULL;

      cloog_upol_set_next (u, n);
      u = n;
      p->next = NULL;
      p = next;
    }

  return res;
}

/**
 * The maximal number of rays allowed to be allocated by PolyLib. In fact since
 * version 5.20, PolyLib automatically tune the number of rays by multiplying
 * by 2 this number each time the maximum is reached. For unknown reasons
 * PolyLib makes a segmentation fault if this number is too small. If this
 * number is too small, performances will be reduced, if it is too high, memory
 * will be saturated. Note that the option "-rays X" set this number to X.
 */
int MAX_RAYS = 50;

/* Unused in this backend.  */

int cloog_domain_allocated = 0;
int cloog_domain_freed = 0;
int cloog_domain_max = 0;

/* The same for Value variables since in GMP mode they have to be freed. */
int cloog_value_allocated = 0;
int cloog_value_freed = 0;
int cloog_value_max = 0;


void
cloog_value_leak_up ()
{
  cloog_value_allocated++;
  if ((cloog_value_allocated - cloog_value_freed) > cloog_value_max)
    cloog_value_max = cloog_value_allocated - cloog_value_freed;
}


void
cloog_value_leak_down ()
{
  cloog_value_freed++;
}

static inline void
cloog_domain_polyhedron_set (CloogDomain * d, ppl_polyhedra_union * p)
{
  d->_polyhedron = p;
}

static inline void
cloog_domain_set_references (CloogDomain * d, int i)
{
  d->_references = i;
}

static CloogDomain *
cloog_new_domain (ppl_polyhedra_union *p)
{
  CloogDomain *domain = (CloogDomain *) malloc (sizeof (CloogDomain));
  domain->_polyhedron = p;
  cloog_domain_set_references (domain, 1);
  return domain;
}

static CloogDomain *
cloog_domain_alloc (Polyhedron *p)
{
  return print_result ("cloog_domain_alloc", cloog_new_domain (p2u (p)));
}

static inline CloogDomain *
cloog_check_domain_id (CloogDomain *dom)
{
  return dom;
}

static inline CloogDomain *
cloog_check_domain (CloogDomain *dom)
{
  if (!dom)
    return dom;

  /* FIXME: Remove this check.  */
  if (cloog_domain_polyhedron (dom)->next)
    {
      fprintf (stderr, "polyhedra of domains should be convex.\n");
      exit (1);
    }

  return dom;
}

/**
 * cloog_domain_matrix2domain function:
 * Given a matrix of constraints (matrix), this function constructs and returns
 * the corresponding domain (i.e. the CloogDomain structure including the
 * polyhedron with its double representation: constraint matrix and the set of
 * rays).
 */
static CloogDomain *
cloog_domain_matrix2domain (CloogMatrix * matrix)
{
  return print_result ("cloog_domain_matrix2domain", cloog_check_domain (cloog_domain_alloc (Constraints2Polyhedron (matrix, MAX_RAYS))));
}

static inline CloogMatrix *
cloog_upol_domain2matrix (ppl_polyhedra_union * upol)
{
  return Polyhedron2Constraints (cloog_upol_polyhedron (upol));
}

/* In the matrix representation an equality has a 0 in the first
   column.  When the value of the first column is 1, the row
   represents an inequality.  */

static inline int
cloog_matrix_row_is_eq_p (CloogMatrix *matrix, int row)
{
  return value_zero_p (matrix->p[row][0]);
}

static ppl_Polyhedron_t
cloog_translate_constraint_matrix (CloogMatrix *matrix)
{
  int i, j;
  ppl_Polyhedron_t res;
  ppl_Constraint_t cstr;
  ppl_Linear_Expression_t expr;
  ppl_Coefficient_t coef;
  ppl_dimension_type dim = matrix->NbColumns - 2;

  ppl_new_Coefficient (&coef);
  ppl_new_NNC_Polyhedron_from_dimension (&res, dim);

  for (i = 0; i < matrix->NbRows; i++)
    {
      ppl_new_Linear_Expression_with_dimension (&expr, dim);

      for (j = 1; j < matrix->NbColumns - 1; j++)
	{
	  ppl_assign_Coefficient_from_mpz_t (coef, matrix->p[i][j]);
	  ppl_Linear_Expression_add_to_coefficient (expr, j - 1, coef);
	}

      ppl_assign_Coefficient_from_mpz_t 
	(coef, matrix->p[i][matrix->NbColumns - 1]);
      ppl_Linear_Expression_add_to_inhomogeneous (expr, coef);

      if (cloog_matrix_row_is_eq_p (matrix, i))
	ppl_new_Constraint (&cstr, expr, PPL_CONSTRAINT_TYPE_EQUAL);
      else
	ppl_new_Constraint (&cstr, expr, PPL_CONSTRAINT_TYPE_GREATER_THAN_OR_EQUAL);

      ppl_Polyhedron_add_constraint (res, cstr);
    }

  if (cloog_check_polyhedral_ops)
    ppl_Polyhedron_OK (res);

  return res;
}

static CloogDomain *
cloog_translate_ppl_polyhedron (ppl_Polyhedron_t pol)
{
  CloogDomain *res;
  CloogMatrix *matrix ;
  ppl_dimension_type dim;
  ppl_const_Constraint_System_t pcs;
  ppl_Constraint_System_const_iterator_t cit, end;
  int row;

  ppl_Polyhedron_constraints (pol, &pcs);
  ppl_new_Constraint_System_const_iterator (&cit);
  ppl_new_Constraint_System_const_iterator (&end);

  for (row = 0, ppl_Constraint_System_begin (pcs, cit), ppl_Constraint_System_end (pcs, end);
       !ppl_Constraint_System_const_iterator_equal_test (cit, end);
       ppl_Constraint_System_const_iterator_increment (cit), row++);

  ppl_Polyhedron_space_dimension (pol, &dim);
  matrix = cloog_matrix_alloc (row, dim + 2);

  for (row = 0, ppl_Constraint_System_begin (pcs, cit), ppl_Constraint_System_end (pcs, end);
       !ppl_Constraint_System_const_iterator_equal_test (cit, end);
       ppl_Constraint_System_const_iterator_increment (cit), row++)
    {
      ppl_const_Constraint_t pc;
      ppl_Coefficient_t coef;
      ppl_dimension_type col;
      Value val;
      int neg;

      value_init (val);
      ppl_new_Coefficient (&coef);
      ppl_Constraint_System_const_iterator_dereference (cit, &pc);

      neg = (ppl_Constraint_type (pc) == PPL_CONSTRAINT_TYPE_LESS_THAN
	     || ppl_Constraint_type (pc) == PPL_CONSTRAINT_TYPE_LESS_THAN_OR_EQUAL) ? 1 : 0;

      for (col = 0; col < dim; col++)
	{
	  ppl_Constraint_coefficient (pc, col, coef);
	  ppl_Coefficient_to_mpz_t (coef, val);

	  if (neg)
	    value_oppose (val, val);

	  value_assign (matrix->p[row][col+1], val);
	}

      ppl_Constraint_inhomogeneous_term (pc, coef);
      ppl_Coefficient_to_mpz_t (coef, val);
      value_assign (matrix->p[row][dim + 1], val);

      switch (ppl_Constraint_type (pc))
	{
	case PPL_CONSTRAINT_TYPE_EQUAL:
	  value_set_si (matrix->p[row][0], 0);
	  break;

	case PPL_CONSTRAINT_TYPE_LESS_THAN:
	case PPL_CONSTRAINT_TYPE_GREATER_THAN:
	  value_decrement (matrix->p[row][dim + 1], matrix->p[row][dim + 1]);
	  value_set_si (matrix->p[row][0], 1);
	  break;

	case PPL_CONSTRAINT_TYPE_LESS_THAN_OR_EQUAL:
	case PPL_CONSTRAINT_TYPE_GREATER_THAN_OR_EQUAL:
	  value_set_si (matrix->p[row][0], 1);
	  break;

	default:
	  fprintf (stderr, "PPL_CONSTRAINT_TYPE_%d not implemented yet\n",
		   ppl_Constraint_type (pc));
	  exit (1);
	}
    }

  res = cloog_domain_matrix2domain (matrix);
  return print_result ("cloog_translate_ppl_polyhedron", cloog_check_domain (res));
}

static inline int
cloog_domain_references (CloogDomain * d)
{
  return d->_references;
}

/**
 * cloog_domain_print function:
 * This function prints the content of a CloogDomain structure (domain) into
 * a file (foo, possibly stdout).
 */
void
cloog_domain_print (FILE * foo, CloogDomain * domain)
{
  ppl_polyhedra_union *upol = cloog_domain_upol (domain);

  while (upol)
    {
      Polyhedron_Print (foo, P_VALUE_FMT, cloog_upol_polyhedron (upol));
      upol = cloog_upol_next (upol);
    }

  fprintf (foo, "Number of active references: %d\n",
	   cloog_domain_references (domain));
}

/**
 * cloog_domain_free function:
 * This function frees the allocated memory for a CloogDomain structure
 * (domain). It decrements the number of active references to this structure,
 * if there are no more references on the structure, it frees it (with the
 * included list of polyhedra).
 */
void
cloog_domain_free (CloogDomain * domain)
{
  if (domain)
    {
      cloog_domain_set_references (domain,
				   cloog_domain_references (domain) - 1);

      if (cloog_domain_references (domain) == 0)
	{

	  ppl_polyhedra_union *upol = cloog_domain_upol (domain);

	  while (upol)
	    {
	      Polyhedron_Free (cloog_upol_polyhedron (upol));
	      upol = cloog_upol_next (upol);
	    }

	  free (domain);
	}
    }
}


/**
 * cloog_domain_copy function:
 * This function returns a copy of a CloogDomain structure (domain). To save
 * memory this is not a memory copy but we increment a counter of active
 * references inside the structure, then return a pointer to that structure.
 */
CloogDomain *
cloog_domain_copy (CloogDomain * domain)
{
  cloog_domain_set_references (domain, cloog_domain_references (domain) + 1);
  return print_result ("cloog_domain_copy", domain);
}


/**
 * cloog_domain_image function:
 * This function returns a CloogDomain structure such that the included
 * polyhedral domain is computed from the former one into another
 * domain according to a given affine mapping function (mapping). 
 */
static CloogDomain *
cloog_domain_image (CloogDomain * domain, CloogMatrix * mapping)
{
  Polyhedron *p = d2p (domain);
  CloogDomain *res =
    cloog_check_domain (cloog_domain_alloc
			(DomainImage (p, mapping, MAX_RAYS)));
  Polyhedron_Free (p);
  return print_result ("cloog_domain_image", res);
}


/**
 * cloog_domain_preimage function:
 * Given a polyhedral domain (polyhedron) inside a CloogDomain structure and a
 * mapping function (mapping), this function returns a new CloogDomain structure
 * with a polyhedral domain which when transformed by mapping function (mapping)
 * gives (polyhedron).
 */
static CloogDomain *
cloog_domain_preimage (CloogDomain * domain, CloogMatrix * mapping)
{
  Polyhedron *p = d2p (domain);
  CloogDomain *res =
    cloog_check_domain (cloog_domain_alloc
			(DomainPreimage (p, mapping, MAX_RAYS)));
  Polyhedron_Free (p);
  return print_result ("cloog_domain_preimage", res);
}


/**
 * cloog_domain_convex function:
 * Given a polyhedral domain (polyhedron), this function concatenates the lists
 * of rays and lines of the two (or more) polyhedra in the domain into one
 * combined list, and  find the set of constraints which tightly bound all of
 * those objects. It returns the corresponding polyhedron.
 */
CloogDomain *
cloog_domain_convex (CloogDomain * domain)
{
  Polyhedron *p = d2p (domain);
  CloogDomain *res =
    cloog_check_domain (cloog_domain_alloc
			(DomainConvex (p, MAX_RAYS)));
  Polyhedron_Free (p);
  return print_result ("cloog_domain_convex", res);
}

static inline unsigned
cloog_upol_nbc (ppl_polyhedra_union * p)
{
  return cloog_upol_polyhedron (p)->NbConstraints;
}

static inline int
cloog_domain_nbconstraints (CloogDomain * domain)
{
  return cloog_domain_polyhedron (domain)->NbConstraints;
}

static inline unsigned
cloog_upol_nbeq (ppl_polyhedra_union * d)
{
  return cloog_upol_polyhedron (d)->NbEq;
}

static inline unsigned
cloog_domain_nbeq (CloogDomain * d)
{
  return cloog_domain_polyhedron (d)->NbEq;
}

static inline unsigned
cloog_upol_dim (ppl_polyhedra_union * p)
{
  return cloog_upol_polyhedron (p)->Dimension;
}

int
cloog_domain_isconvex (CloogDomain * domain)
{
  if (cloog_domain_polyhedron (domain))
    return !cloog_upol_next (cloog_domain_upol (domain));

  return 1;
}

unsigned
cloog_domain_dim (CloogDomain * d)
{
  return cloog_domain_polyhedron (d)->Dimension;
}

static CloogDomain *
cloog_check_domains (CloogDomain *ppl, CloogDomain *polylib)
{
  /* Cannot use cloog_domain_lazy_equal (polylib, ppl) here as this
     function is too dumb: it does not detect permutations of
     constraints.  */
  if (!cloog_domain_isempty (cloog_domain_difference (ppl, polylib))
      || !cloog_domain_isempty (cloog_domain_difference (polylib, ppl)))
    {
      fprintf (stderr, "different domains ( \n ppl (\n");
      cloog_domain_print (stderr, ppl);
      fprintf (stderr, ") \n polylib (\n");
      cloog_domain_print (stderr, polylib);
      fprintf (stderr, "))\n");
      exit (1);
    }

  if (cloog_return_ppl_result)
    return ppl;
  else
    return polylib;
}

/**
 * cloog_domain_simple_convex:
 * Given a list (union) of polyhedra, this function returns a "simple"
 * convex hull of this union.  In particular, the constraints of the
 * the returned polyhedron consist of (parametric) lower and upper
 * bounds on individual variables and constraints that appear in the
 * original polyhedra.
 *
 * nb_par is the number of parameters of the domain.
 */
CloogDomain *
cloog_domain_simple_convex (CloogDomain * domain, int nb_par)
{
  fprintf (stderr, "cloog_domain_simple_convex is not implemented yet.\n");
  exit (1);
}


/**
 * cloog_domain_simplify function:
 * Given two polyhedral domains (pol1) and (pol2) inside two CloogDomain
 * structures, this function finds the largest domain set (or the smallest list
 * of non-redundant constraints), that when intersected with polyhedral
 * domain (pol2) equals (Pol1)intersect(Pol2). The output is a new CloogDomain
 * structure with a polyhedral domain with the "redundant" constraints removed.
 * NB: this function do not work as expected with unions of polyhedra...
 */
CloogDomain *
cloog_domain_simplify (CloogDomain * dom1, CloogDomain * dom2)
{
  CloogMatrix *M, *M2;
  CloogDomain *dom;
  Polyhedron *d2, *P = d2p (dom1);
  int nbc = cloog_domain_nbconstraints (dom1);

  /* DomainSimplify doesn't remove all redundant equalities,
   * so we remove them here first in case both dom1 and dom2
   * are single polyhedra (i.e., not unions of polyhedra).
   */
  if (cloog_domain_isconvex (dom1) && cloog_domain_isconvex (dom2)
      && cloog_domain_nbeq (dom1) && cloog_domain_nbeq (dom2))
    {
      int i, row;
      int rows = cloog_domain_nbeq (dom1) + cloog_domain_nbeq (dom2);
      int cols = cloog_domain_dim (dom1) + 2;
      int rank;
      M = cloog_matrix_alloc (rows, cols);
      M2 = cloog_matrix_alloc (nbc, cols);
      Vector_Copy (cloog_domain_polyhedron (dom2)->Constraint[0],
		   M->p[0], cloog_domain_nbeq (dom2) * cols);
      rank = cloog_domain_nbeq (dom2);
      row = 0;
      for (i = 0; i < cloog_domain_nbeq (dom1); ++i)
	{
	  Vector_Copy (P->Constraint[i], M->p[rank], cols);
	  if (Gauss (M, rank + 1, cols - 1) > rank)
	    {
	      Vector_Copy (P->Constraint[i], M2->p[row++], cols);
	      rank++;
	    }
	}
      if (row < cloog_domain_nbeq (dom1))
	{
	  Vector_Copy (P->Constraint[cloog_domain_nbeq (dom1)],
		       M2->p[row], (nbc - cloog_domain_nbeq (dom1)) * cols);
	  P = Constraints2Polyhedron (M2, MAX_RAYS);
	}
      cloog_matrix_free (M2);
      cloog_matrix_free (M);
    }
  d2 = d2p (dom2);
  dom = cloog_domain_alloc (DomainSimplify (P, d2, MAX_RAYS));
  Polyhedron_Free (d2);
  Polyhedron_Free (P);
  return print_result ("cloog_domain_simplify", cloog_check_domain (dom));
}

static ppl_polyhedra_union *
cloog_upol_copy (ppl_polyhedra_union *p)
{
  ppl_polyhedra_union *res = cloog_new_upol (Polyhedron_Copy (cloog_upol_polyhedron (p)));
  ppl_polyhedra_union *upol = res;

  while (cloog_upol_next (p))
    {
      cloog_upol_set_next (upol, cloog_new_upol (Polyhedron_Copy (cloog_upol_polyhedron (p))));
      upol = cloog_upol_next (upol);
      p = cloog_upol_next (p);
    }

  return res;
}

/* Adds to D1 the union of polyhedra from D2, with no checks of
   redundancies between polyhedra.  */

CloogDomain *
cloog_domain_add_domain (CloogDomain *d1, CloogDomain *d2)
{
  ppl_polyhedra_union *upol;

  if (!d1)
    return d2;

  if (!d2)
    return d1;

  upol = cloog_domain_upol (d1);
  while (cloog_upol_next (upol))
    upol = cloog_upol_next (upol);

  cloog_upol_set_next (upol, cloog_domain_upol (d2));
  return d1;
}

/**
 * cloog_domain_union function:
 * This function returns a new CloogDomain structure including a polyhedral
 * domain which is the union of two polyhedral domains (pol1) U (pol2) inside
 * two CloogDomain structures.
 */
CloogDomain *
cloog_domain_union (CloogDomain * dom1, CloogDomain * dom2)
{
  CloogDomain *res;
  ppl_polyhedra_union *head1, *head2, *tail1, *tail2;
  ppl_polyhedra_union *d1, *d2;

  if (!dom1)
    return dom2;

  if (!dom2)
    return dom1;

  if (cloog_domain_dim (dom1) != cloog_domain_dim (dom2))
    {
      fprintf (stderr, "cloog_domain_union should not be called on domains of different dimensions.\n");
      exit (1);
    }

  head1 = NULL;
  tail1 = NULL;
  for (d1 = cloog_domain_upol (dom1); d1; d1 = cloog_upol_next (d1))
    {
      int found = 0;
      ppl_Polyhedron_t ppl1 = cloog_translate_constraint_matrix (cloog_upol_domain2matrix (d1));

      for (d2 = cloog_domain_upol (dom2); d2; d2 = cloog_upol_next (d2))
	{
	  ppl_Polyhedron_t ppl2 = cloog_translate_constraint_matrix (cloog_upol_domain2matrix (d2));

	  if (ppl_Polyhedron_contains_Polyhedron (ppl2, ppl1))
	    {
	      found = 1;
	      break;
	    }
	}

      if (!found)
	{
	  if (!tail1)
	    {
	      head1 = cloog_upol_copy (d1);
	      tail1 = head1;
	    }
	  else
	    {
	      cloog_upol_set_next (tail1, cloog_upol_copy (d1));
	      tail1 = cloog_upol_next (tail1);
	    }
	}
    }

  head2 = NULL;
  tail2 = NULL;
  for (d2 = cloog_domain_upol (dom2); d2; d2 = cloog_upol_next (d2))
    {
      int found = 0;
      ppl_Polyhedron_t ppl2 = cloog_translate_constraint_matrix (cloog_upol_domain2matrix (d2));

      for (d1 = head1; d1; d1 = cloog_upol_next (d1))
	{
	  ppl_Polyhedron_t ppl1 = cloog_translate_constraint_matrix (cloog_upol_domain2matrix (d1));

	  if (ppl_Polyhedron_contains_Polyhedron (ppl1, ppl2))
	    {
	      found = 1;
	      break;
	    }
	}

      if (!found)
	{
	  if (!tail2)
	    {
	      head2 = cloog_upol_copy (d2);
	      tail2 = head2;
	    }
	  else
	    {
	      cloog_upol_set_next (tail2, cloog_upol_copy (d2));
	      tail2 = cloog_upol_next (tail2);
	    }
	}
    }

  if (!head1)
    res = cloog_new_domain (head2);
  else
    {
      cloog_upol_set_next (tail1, head2);
      res = cloog_new_domain (head1);
    }

  if (cloog_check_polyhedral_ops)
    {
      Polyhedron *p1 = d2p (dom1);
      Polyhedron *p2 = d2p (dom2);

      cloog_check_domains (res, cloog_domain_alloc (DomainUnion (p1, p2, MAX_RAYS)));

      Polyhedron_Free (p1);
      Polyhedron_Free (p2);
    }

  return print_result ("cloog_domain_union", cloog_check_domain (res));
}

/**
 * cloog_domain_intersection function:
 * This function returns a new CloogDomain structure including a polyhedral
 * domain which is the intersection of two polyhedral domains (pol1)inter(pol2)
 * inside two CloogDomain structures.
 */
CloogDomain *
cloog_domain_intersection (CloogDomain * dom1, CloogDomain * dom2)
{
  CloogDomain *res;
  ppl_polyhedra_union *p1, *p2;
  ppl_Polyhedron_t ppl1, ppl2;

  res = cloog_domain_empty (cloog_domain_dim (dom1));

  for (p1 = cloog_domain_upol (dom1); p1; p1 = cloog_upol_next (p1))
    {
      ppl1 = cloog_translate_constraint_matrix (Polyhedron2Constraints (cloog_upol_polyhedron (p1)));

      for (p2 = cloog_domain_upol (dom2); p2; p2 = cloog_upol_next (p2))
	{
	  ppl2 = cloog_translate_constraint_matrix (Polyhedron2Constraints (cloog_upol_polyhedron (p2)));
	  ppl_Polyhedron_intersection_assign (ppl2, ppl1);

	  res = cloog_domain_union (res, cloog_translate_ppl_polyhedron (ppl2));
	}
    }

  if (cloog_check_polyhedral_ops)
    {
      Polyhedron *a1 = d2p (dom1);
      Polyhedron *a2 = d2p (dom2);

      res = cloog_check_domains (res, cloog_domain_alloc (DomainIntersection (a1, a2, MAX_RAYS)));

      Polyhedron_Free (a1);
      Polyhedron_Free (a2);
    }

  return print_result ("cloog_domain_intersection", res);
}


/**
 * cloog_domain_difference function:
 * This function returns a new CloogDomain structure including a polyhedral
 * domain which is the difference of two polyhedral domains domain \ minus
 * inside two CloogDomain structures.
 * - November 8th 2001: first version.
 */
CloogDomain *
cloog_domain_difference (CloogDomain * domain, CloogDomain * minus)
{
  if (cloog_domain_isempty (minus))
    return print_result ("cloog_domain_difference", cloog_domain_copy (domain));
  else
    {
      Polyhedron *p1 = d2p (domain);
      Polyhedron *p2 = d2p (minus);
      CloogDomain *res = cloog_domain_alloc (DomainDifference (p1, p2, MAX_RAYS));
      Polyhedron_Free (p1);
      Polyhedron_Free (p2);
      return print_result ("cloog_domain_difference", res);
    }
}


/**
 * cloog_domain_addconstraints function :
 * This function adds source's polyhedron constraints to target polyhedron: for
 * each element of the polyhedron inside 'target' (i.e. element of the union
 * of polyhedra) it adds the constraints of the corresponding element in
 * 'source'.
 * - August 10th 2002: first version.
 * Nota bene for future : it is possible that source and target don't have the
 * same number of elements (try iftest2 without non-shared constraint
 * elimination in cloog_loop_separate !). This function is yet another part
 * of the DomainSimplify patching problem...
 */
CloogDomain *
cloog_domain_addconstraints (domain_source, domain_target)
     CloogDomain *domain_source, *domain_target;
{
  unsigned nb_constraint;
  Value *constraints;
  ppl_polyhedra_union *source, *target, *new, *next, *last;

  source = cloog_domain_upol (domain_source);
  target = cloog_domain_upol (domain_target);

  constraints = cloog_upol_polyhedron (source)->p_Init;
  nb_constraint = cloog_upol_nbc (source);
  last = new = cloog_new_upol (AddConstraints (constraints, nb_constraint,
					       u2p (target), MAX_RAYS));
  source = cloog_upol_next (source);
  next = cloog_upol_next (target);

  while (next)
    {				/* BUG !!! This is actually a bug. I don't know yet how to cleanly avoid
				 * the situation where source and target do not have the same number of
				 * elements. So this 'if' is an awful trick, waiting for better.
				 */
      if (source)
	{
	  constraints = cloog_upol_polyhedron (source)->p_Init;
	  nb_constraint = cloog_upol_nbc (source);
	  source = cloog_upol_next (source);
	}
      cloog_upol_set_next 
	(last, cloog_new_upol (AddConstraints (constraints, nb_constraint,
					       u2p (next), MAX_RAYS)));
      last = cloog_upol_next (last);
      next = cloog_upol_next (next);
    }

  return print_result ("cloog_domain_addconstraints", cloog_check_domain (cloog_new_domain (new)));
}


/**
 * cloog_domain_sort function:
 * This function topologically sorts (nb_pols) polyhedra. Here (pols) is a an
 * array of pointers to polyhedra, (nb_pols) is the number of polyhedra,
 * (level) is the level to consider for partial ordering (nb_par) is the
 * parameter space dimension, (permut) if not NULL, is an array of (nb_pols)
 * integers that contains a permutation specification after call in order to
 * apply the topological sorting. 
 */
void
cloog_domain_sort (doms, nb_pols, level, nb_par, permut)
     CloogDomain **doms;
     unsigned nb_pols, level, nb_par;
     int *permut;
{
  int *time, i;
  Polyhedron **pols = (Polyhedron **) malloc (nb_pols * sizeof (Polyhedron *));

  for (i = 0; i < nb_pols; i++)
    pols[i] = cloog_domain_polyhedron (doms[i]);

  /* time is an array of (nb_pols) integers to store logical time values. We
   * do not use it, but it is compulsory for PolyhedronTSort.
   */
  time = (int *) malloc (nb_pols * sizeof (int));

  /* PolyhedronTSort will fill up permut (and time). */
  PolyhedronTSort (pols, nb_pols, level, nb_par, time, permut, MAX_RAYS);

  free (pols);
  free (time);
}


/**
 * cloog_domain_empty function:
 * This function allocates the memory space for a CloogDomain structure and
 * sets its polyhedron to an empty polyhedron with 'dimension' dimensions.
 * Then it returns a pointer to the allocated space.
 * - June 10th 2005: first version.
 */
CloogDomain *
cloog_domain_empty (int dimension)
{
  return (cloog_domain_alloc (Empty_Polyhedron (dimension)));
}


/******************************************************************************
 *                          Structure display function                        *
 ******************************************************************************/


/**
 * cloog_domain_print_structure :
 * this function is a more human-friendly way to display the CloogDomain data
 * structure, it only shows the constraint system and includes an indentation
 * level (level) in order to work with others print_structure functions.
 * Written by Olivier Chorier, Luc Marchaud, Pierre Martin and Romain Tartiere.
 * - April 24th 2005: Initial version.
 * - May   26th 2005: Memory leak hunt.
 * - June  16th 2005: (Ced) Integration in domain.c.
 */
void
cloog_domain_print_structure (FILE * file, CloogDomain * domain, int level)
{
  int i;
  CloogMatrix *matrix;

  /* Go to the right level. */
  for (i = 0; i < level; i++)
    fprintf (file, "|\t");

  if (domain != NULL)
    {
      fprintf (file, "+-- CloogDomain\n");

      /* Print the matrix. */
      matrix = cloog_upol_domain2matrix (cloog_domain_upol (domain));
      cloog_matrix_print_structure (file, matrix, level);
      cloog_matrix_free (matrix);

      /* A blank line. */
      for (i = 0; i < level + 1; i++)
	fprintf (file, "|\t");
      fprintf (file, "\n");
    }
  else
    fprintf (file, "+-- Null CloogDomain\n");

}


/**
 * cloog_domain_list_print function:
 * This function prints the content of a CloogDomainList structure into a
 * file (foo, possibly stdout).
 * - November 6th 2001: first version.
 */
void
cloog_domain_list_print (FILE * foo, CloogDomainList * list)
{
  while (list != NULL)
    {
      cloog_domain_print (foo, cloog_domain (list));
      list = cloog_next_domain (list);
    }
}


/******************************************************************************
 *                         Memory deallocation function                       *
 ******************************************************************************/


/**
 * cloog_domain_list_free function:
 * This function frees the allocated memory for a CloogDomainList structure.
 * - November 6th 2001: first version.
 */
void
cloog_domain_list_free (CloogDomainList * list)
{
  CloogDomainList *temp;

  while (list != NULL)
    {
      temp = cloog_next_domain (list);
      cloog_domain_free (cloog_domain (list));
      free (list);
      list = temp;
    }
}


/******************************************************************************
 *                               Reading function                             *
 ******************************************************************************/


/**
 * cloog_domain_read function:
 * Adaptation from the PolyLib. This function reads a matrix into a file (foo,
 * posibly stdin) and returns a pointer to a polyhedron containing the read
 * information. 
 * - October 18th 2001: first version.
 */
CloogDomain *
cloog_domain_read (FILE * foo)
{
  CloogMatrix *matrix;
  CloogDomain *domain;

  matrix = cloog_matrix_read (foo);
  domain = cloog_domain_matrix2domain (matrix);
  cloog_matrix_free (matrix);

  return print_result ("cloog_domain_read", domain);
}


/**
 * cloog_domain_union_read function:
 * This function reads a union of polyhedra into a file (foo, posibly stdin) and
 * returns a pointer to a Polyhedron containing the read information. 
 * - September 9th 2002: first version.
 * - October  29th 2005: (debug) removal of a leak counting "correction" that
 *                       was just false since ages.
 */
CloogDomain *
cloog_domain_union_read (FILE * foo)
{
  int i, nb_components;
  char s[MAX_STRING];
  CloogDomain *domain, *temp, *old;

  /* domain reading: nb_components (constraint matrices). */
  while (fgets (s, MAX_STRING, foo) == 0);
  while ((*s == '#' || *s == '\n') || (sscanf (s, " %d", &nb_components) < 1))
    fgets (s, MAX_STRING, foo);

  if (nb_components > 0)
    {				/* 1. first part of the polyhedra union, */
      domain = cloog_domain_read (foo);
      /* 2. and the nexts. */
      for (i = 1; i < nb_components; i++)
	{			/* Leak counting is OK since next allocated domain is freed here. */
	  temp = cloog_domain_read (foo);
	  old = domain;
	  domain = cloog_domain_union (temp, old);
	  cloog_domain_free (temp);
	  cloog_domain_free (old);
	}
      return print_result ("cloog_domain_union_read", cloog_check_domain (domain));
    }
  else
    return NULL;
}


/**
 * cloog_domain_list_read function:
 * This function reads a list of polyhedra into a file (foo, posibly stdin) and
 * returns a pointer to a CloogDomainList containing the read information. 
 * - November 6th 2001: first version.
 */
CloogDomainList *
cloog_domain_list_read (FILE * foo)
{
  int i, nb_pols;
  char s[MAX_STRING];
  CloogDomainList *list, *now, *next;


  /* We read first the number of polyhedra in the list. */
  while (fgets (s, MAX_STRING, foo) == 0);
  while ((*s == '#' || *s == '\n') || (sscanf (s, " %d", &nb_pols) < 1))
    fgets (s, MAX_STRING, foo);

  /* Then we read the polyhedra. */
  list = NULL;
  if (nb_pols > 0)
    {
      list = (CloogDomainList *) malloc (sizeof (CloogDomainList));
      cloog_set_domain (list, cloog_domain_read (foo));
      cloog_set_next_domain (list, NULL);
      now = list;
      for (i = 1; i < nb_pols; i++)
	{
	  next = (CloogDomainList *) malloc (sizeof (CloogDomainList));
	  cloog_set_domain (next, cloog_domain_read (foo));
	  cloog_set_next_domain (next, NULL);
	  cloog_set_next_domain (now, next);
	  now = cloog_next_domain (now);
	}
    }
  return list;
}


/******************************************************************************
 *                            Processing functions                            *
 ******************************************************************************/

/**
 * cloog_domain_isempty function:
 * This function returns 1 if the polyhedron given as input is empty, 0
 * otherwise.
 * - October 28th 2001: first version.
 */
int
cloog_domain_isempty (CloogDomain * domain)
{
  if (cloog_domain_polyhedron (domain) == NULL)
    return 1;

  if (cloog_upol_next (cloog_domain_upol (domain)))
    return 0;

  return ((cloog_domain_dim (domain) < cloog_domain_nbeq (domain)) ? 1 : 0);
}

/**
 * cloog_domain_project function:
 * From Quillere's LoopGen 0.4. This function returns the projection of
 * (domain) on the (level) first dimensions (i.e. outer loops). It returns a
 * pointer to the projected Polyhedron.
 * - nb_par is the number of parameters.
 **
 * - October 27th 2001: first version.
 * - June 21rd 2005: Adaptation for GMP (based on S. Verdoolaege's version of
 *                   CLooG 0.12.1).
 */
CloogDomain *
cloog_domain_project_1 (CloogDomain * domain, int level, int nb_par)
{
  int row, column, nb_rows, nb_columns, difference;
  CloogDomain *projected_domain;
  CloogMatrix *matrix;

  nb_rows = level + nb_par + 1;
  nb_columns = cloog_domain_dim (domain) + 1;
  difference = nb_columns - nb_rows;

  if (difference == 0)
    return print_result ("cloog_domain_project", cloog_domain_copy (domain));

  matrix = cloog_matrix_alloc (nb_rows, nb_columns);

  for (row = 0; row < level; row++)
    for (column = 0; column < nb_columns; column++)
      value_set_si (matrix->p[row][column], (row == column ? 1 : 0));

  for (; row < nb_rows; row++)
    for (column = 0; column < nb_columns; column++)
      value_set_si (matrix->p[row][column],
		    (row + difference == column ? 1 : 0));

  projected_domain = cloog_domain_image (domain, matrix);
  cloog_matrix_free (matrix);

  return print_result ("cloog_domain_project_1", cloog_check_domain (projected_domain));
}

CloogDomain *
cloog_domain_project (CloogDomain * domain, int level, int nb_par)
{
  return print_result ("cloog_domain_project",
		       cloog_domain_project_1 (domain, level, nb_par));
}

CloogDomain *
cloog_domain_project_ported (CloogDomain * domain, int level, int nb_par)
{
  CloogDomain *res = NULL;
  ppl_polyhedra_union *upol = cloog_domain_upol (domain);
  int i, diff = cloog_domain_dim (domain) - level - nb_par;
  size_t n;
  ppl_dimension_type *to_remove;

  if (diff == 0)
    return print_result ("cloog_domain_project", domain);

  n = diff;
  to_remove = (ppl_dimension_type *) malloc (n * sizeof (ppl_dimension_type));

  for (i = 0; i < n; i++)
    to_remove[i] = i + level;

  while (upol)
    {
      ppl_Polyhedron_t ppl = cloog_translate_constraint_matrix (cloog_upol_domain2matrix (upol));

      ppl_Polyhedron_remove_space_dimensions (ppl, to_remove, n);
      res = cloog_domain_add_domain (res, cloog_translate_ppl_polyhedron (ppl));
      upol = cloog_upol_next (upol);
    }

  if (cloog_check_polyhedral_ops)
    return print_result ("cloog_domain_project", 
			 cloog_check_domains
			 (res, cloog_domain_project_1 (domain, level, nb_par)));

  return print_result ("cloog_domain_project", res);
}

/**
 * cloog_domain_extend function:
 * From Quillere's LoopGen 0.4. This function returns the (domain) given as
 * input with (dim)+(nb_par) dimensions. The new dimensions are added before
 * the (nb_par) parameters. This function does not free (domain), and returns
 * a new polyhedron.
 * - nb_par is the number of parameters.
 **
 * - October 27th 2001: first version.
 * - June 21rd 2005: Adaptation for GMP (based on S. Verdoolaege's version of
 *                 CLooG 0.12.1).
 */
CloogDomain *
cloog_domain_extend (CloogDomain * domain, int dim, int nb_par)
{
  int row, column, nb_rows, nb_columns, difference;
  CloogDomain *extended_domain;
  CloogMatrix *matrix;

  nb_rows = 1 + cloog_domain_dim (domain);
  nb_columns = dim + nb_par + 1;
  difference = nb_columns - nb_rows;

  if (difference == 0)
    return print_result ("cloog_domain_extend", cloog_domain_copy (domain));

  matrix = cloog_matrix_alloc (nb_rows, nb_columns);

  for (row = 0; row < cloog_domain_dim (domain) - nb_par; row++)
    for (column = 0; column < nb_columns; column++)
      value_set_si (matrix->p[row][column], (row == column ? 1 : 0));

  for (; row <= cloog_domain_dim (domain); row++)
    for (column = 0; column < nb_columns; column++)
      value_set_si (matrix->p[row][column],
		    (row + difference == column ? 1 : 0));

  extended_domain = cloog_domain_preimage (domain, matrix);
  cloog_matrix_free (matrix);

  return print_result ("cloog_domain_extend", cloog_check_domain (extended_domain));
}


/**
 * cloog_domain_never_integral function:
 * For us, an equality like 3*i -4 = 0 is always false since 4%3 != 0. This
 * function returns a boolean set to 1 if there is this kind of 'never true'
 * constraint inside a polyhedron, 0 otherwise.
 * - domain is the polyhedron to check,
 **
 * - November 28th 2001: first version. 
 * - June 26th 2003: for iterators, more 'never true' constraints are found
 *                   (compare cholesky2 and vivien with a previous version),
 *                   checking for the parameters created (compare using vivien).
 * - June 28th 2003: Previously in loop.c and called
 *                   cloog_loop_simplify_nevertrue, now here !
 * - June 21rd 2005: Adaptation for GMP (based on S. Verdoolaege's version of
 *                   CLooG 0.12.1).
 * - October 14th 2005: Complete rewriting, not faster but code quite shorter.
 */
int
cloog_domain_never_integral (CloogDomain * domain)
{
  int i, dimension, nbc;
  Value gcd, modulo;
  Polyhedron *polyhedron;

  if ((domain == NULL) || (cloog_domain_polyhedron (domain) == NULL))
    return 1;

  value_init_c (gcd);
  value_init_c (modulo);
  polyhedron = d2p (domain);
  dimension = cloog_domain_dim (domain) + 2;
  nbc = cloog_domain_nbconstraints (domain);

  /* For each constraint... */
  for (i = 0; i < nbc; i++)
    {				/* If we have an equality and the scalar part is not zero... */
      if (value_zero_p (polyhedron->Constraint[i][0]) &&
	  value_notzero_p (polyhedron->Constraint[i][dimension - 1]))
	{			/* Then we check whether the scalar can be divided by the gcd of the
				 * unknown vector (including iterators and parameters) or not. If not,
				 * there is no integer point in the polyhedron and we return 1.
				 */
	  Vector_Gcd (&(polyhedron->Constraint[i][1]), dimension - 2, &gcd);
	  value_modulus (modulo,
			 polyhedron->Constraint[i][dimension - 1],
			 gcd);

	  if (value_notzero_p (modulo))
	    {
	      value_clear_c (gcd);
	      value_clear_c (modulo);
	      Polyhedron_Free (polyhedron);
	      return 1;
	    }
	}
    }

  value_clear_c (gcd);
  value_clear_c (modulo);
  Polyhedron_Free (polyhedron);
  return (0);
}


/**
 * cloog_domain_stride function:
 * This function finds the stride imposed to unknown with the column number
 * 'strided_level' in order to be integral. For instance, if we have a
 * constraint like -i - 2j + 2k = 0, and we consider k, then k can be integral
 * only if (i + 2j)%2 = 0. Then only if i%2 = 0. Then k imposes a stride 2 to
 * the unknown i. The function returns the imposed stride in a parameter field.
 * - domain is the set of constraint we have to consider,
 * - strided_level is the column number of the unknown for which a stride have
 *   to be found,
 * - looking_level is the column number of the unknown that impose a stride to
 *   the first unknown.
 * - stride is the stride that is returned back as a function parameter. 
 * - offset is the value of the constant c if the condition is of the shape
 *   (i + c)%s = 0, s being the stride.
 **
 * - June 28th 2003: first version.
 * - July 14th 2003: can now look for multiple striding constraints and returns
 *                   the GCD of the strides and the common offset.
 * - June 21rd 2005: Adaptation for GMP (based on S. Verdoolaege's version of
 *                   CLooG 0.12.1).
 */
void
cloog_domain_stride (domain, strided_level, nb_par, stride, offset)
     CloogDomain *domain;
     int strided_level, nb_par;
     Value *stride, *offset;
{
  int i;
  int n_col, n_row, rank;
  CloogMatrix *M;
  Matrix *U;
  Vector *V;
  Polyhedron *polyhedron = d2p (domain);
  int dimension = cloog_domain_dim (domain);
  int nbeq = cloog_domain_nbeq (domain);

  /* Look at all equalities involving strided_level and the inner
   * iterators.  We can ignore the outer iterators and the parameters
   * here because the lower bound on strided_level is assumed to
   * be a constant.
   */
  n_col = (1 + dimension - nb_par) - strided_level;
  for (i = 0, n_row = 0; i < nbeq; i++)
    if (First_Non_Zero
	(polyhedron->Constraint[i] + strided_level, n_col) != -1)
      ++n_row;

  M = cloog_matrix_alloc (n_row + 1, n_col + 1);
  for (i = 0, n_row = 0; i < nbeq; i++)
    {
      if (First_Non_Zero
	  (polyhedron->Constraint[i] + strided_level, n_col) == -1)
	continue;
      Vector_Copy (polyhedron->Constraint[i] + strided_level,
		   M->p[n_row], n_col);
      value_assign (M->p[n_row][n_col],
		    polyhedron->Constraint[i][1 + dimension]);
      ++n_row;
    }
  value_set_si (M->p[n_row][n_col], 1);

  /* Then look at the general solution to the above equalities. */
  rank = SolveDiophantine (M, &U, &V);
  cloog_matrix_free (M);

  if (rank == -1)
    {
      /* There is no solution, so the body of this loop will
       * never execute.  We just leave the constraints alone here so
       * that they will ensure the body will not be executed.
       * We should probably propagate this information up so that
       * the loop can be removed entirely.
       */
      value_set_si (*offset, 0);
      value_set_si (*stride, 1);
    }
  else
    {
      /* Compute the gcd of the coefficients defining strided_level. */
      Vector_Gcd (U->p[0], U->NbColumns, stride);
      value_oppose (*offset, V->p[0]);
      value_pmodulus (*offset, *offset, *stride);
    }
  Matrix_Free (U);
  Vector_Free (V);
  Polyhedron_Free (polyhedron);
  return;
}


/**
 * cloog_domain_integral_lowerbound function:
 * This function returns 1 if the lower bound of an iterator (such as its
 * column rank in the constraint set 'domain' is 'level') is integral,
 * 0 otherwise. If the lower bound is actually integral, the function fills
 * the 'lower' field with the lower bound value.
 * - June 29th 2003: first version.
 * - June 21rd 2005: Adaptation for GMP (based on S. Verdoolaege's version of
 *                   CLooG 0.12.1).
 */
int
cloog_domain_integral_lowerbound (domain, level, lower)
     CloogDomain *domain;
     int level;
     Value *lower;
{
  int i, first_lower = 1, dimension, lower_constraint = -1, nbc;
  Value iterator, constant, tmp;
  Polyhedron *polyhedron;

  polyhedron = d2p (domain);
  dimension = cloog_domain_dim (domain);
  nbc = cloog_domain_nbconstraints (domain);

  /* We want one and only one lower bound (e.g. no equality, no maximum
   * calculation...).
   */
  for (i = 0; i < nbc; i++)
    if (value_zero_p (polyhedron->Constraint[i][0])
	&& value_notzero_p (polyhedron->Constraint[i][level]))
      {
	Polyhedron_Free (polyhedron);
	return 0;
      }

  for (i = 0; i < nbc; i++)
    if (value_pos_p (polyhedron->Constraint[i][level]))
      {
	if (first_lower)
	  {
	    first_lower = 0;
	    lower_constraint = i;
	  }
	else
	  {
	    Polyhedron_Free (polyhedron);
	    return 0;
	  }
      }

  if (first_lower)
    {
      Polyhedron_Free (polyhedron);
      return 0;
    }

  /* We want an integral lower bound: no other non-zero entry except the
   * iterator coefficient and the constant.
   */
  for (i = 1; i < level; i++)
    if (value_notzero_p
	(polyhedron->Constraint[lower_constraint][i]))
      {
	Polyhedron_Free (polyhedron);
	return 0;
      }

  for (i = level + 1; i <= dimension; i++)
    if (value_notzero_p
	(polyhedron->Constraint[lower_constraint][i]))
      {
	Polyhedron_Free (polyhedron);
	return 0;
      }

  value_init_c (iterator);
  value_init_c (constant);
  value_init_c (tmp);

  /* If all is passed, then find the lower bound and return 1. */
  value_assign (iterator,
		polyhedron->Constraint[lower_constraint][level]);
  value_oppose (constant,
		polyhedron->Constraint[lower_constraint][dimension + 1]);

  value_modulus (tmp, constant, iterator);
  value_division (*lower, constant, iterator);

  if (!(value_zero_p (tmp) || value_neg_p (constant)))
    value_increment (*lower, *lower);

  value_clear_c (iterator);
  value_clear_c (constant);
  value_clear_c (tmp);
  Polyhedron_Free (polyhedron);
  return 1;
}


/**
 * cloog_domain_lowerbound_update function:
 * This function updates the integral lower bound of an iterator (such as its
 * column rank in the constraint set 'domain' is 'level') into  'lower'.
 * - Jun  29th 2003: first version.
 * - June 21rd 2005: Adaptation for GMP (based on S. Verdoolaege's version of
 *                   CLooG 0.12.1).
 */
void
cloog_domain_lowerbound_update (domain, level, lower)
     CloogDomain *domain;
     int level;
     Value lower;
{
  int i;
  int nbc = cloog_domain_nbconstraints (domain);
  int dim = cloog_domain_dim (domain);
  Polyhedron *polyhedron = cloog_domain_polyhedron (domain);

  /* There is only one lower bound, the first one is the good one. */
  for (i = 0; i < nbc; i++)
    if (value_pos_p (polyhedron->Constraint[i][level]))
      {
	value_set_si (polyhedron->Constraint[i][level], 1);
	value_oppose (polyhedron->Constraint[i][dim + 1], lower);
	break;
      }
}


/**
 * cloog_domain_lazy_equal function:
 * This function returns 1 if the domains given as input are the same, 0 if it
 * is unable to decide. This function makes an entry-to-entry comparison between
 * the constraint systems, if all the entries are the same, the domains are
 * obviously the same and it returns 1, at the first difference, it returns 0.
 * This is a very fast way to verify this property. It has been shown (with the
 * CLooG benchmarks) that operations on equal domains are 17% of all the
 * polyhedral computations. For 75% of the actually identical domains, this
 * function answer that they are the same and allow to give immediately the
 * trivial solution instead of calling the heavy general functions of PolyLib.
 * - August 22th 2003: first version.
 * - June 21rd 2005: Adaptation for GMP (based on S. Verdoolaege's version of
 *                   CLooG 0.12.1).
 */
int
cloog_domain_lazy_equal (CloogDomain * d1, CloogDomain * d2)
{
  int i, nb_elements;
  ppl_polyhedra_union *u1 = cloog_domain_upol (d1);
  ppl_polyhedra_union *u2 = cloog_domain_upol (d2);

  while (u1 && u2)
    {
      if ((cloog_upol_nbc (u1) != cloog_upol_nbc (u2)) ||
	  (cloog_upol_dim (u1) != cloog_upol_dim (u2)))
	return 0;

      nb_elements =
	cloog_upol_nbc (u1) * (cloog_upol_dim (u1) + 2);

      for (i = 0; i < nb_elements; i++)
	if (value_ne (cloog_upol_polyhedron (u1)->p_Init[i],
		      cloog_upol_polyhedron (u2)->p_Init[i]))
	  return 0;

      u1 = cloog_upol_next (u1);
      u2 = cloog_upol_next (u2);
    }

  if (u1 || u2)
    return 0;

  return 1;
}


/**
 * cloog_domain_lazy_block function:
 * This function returns 1 if the two domains d1 and d2 given as input are the
 * same (possibly except for a dimension equal to a constant where we accept
 * a difference of 1) AND if we are sure that there are no other domain in
 * the code generation problem that may put integral points between those of
 * d1 and d2 (0 otherwise). In fact this function answers the question "can I
 * safely consider the two domains as only one with two statements (a block) ?".
 * This function is lazy: it asks for very standard scattering representation
 * (only one constraint per dimension which is an equality, and the constraints
 * are ordered per dimension depth: the left hand side of the constraint matrix
 * is the identity) and will answer NO at the very first problem.
 * - d1 and d2 are the two domains to check for blocking,
 * - scattering is the linked list of all domains,
 * - scattdims is the total number of scattering dimentions.
 **
 * - April   30th 2005: beginning
 * - June     9th 2005: first working version.
 * - June    10th 2005: debugging.
 * - June    21rd 2005: Adaptation for GMP.
 * - October 16th 2005: (debug) some false blocks have been removed.
 */
int
cloog_domain_lazy_block (d1, d2, scattering, scattdims)
     CloogDomain *d1, *d2;
     CloogDomainList *scattering;
     int scattdims;
{
  int i, j, difference = 0, different_constraint = 0, nbc;
  int dim1, dim2;
  Value date1, date2, date3, temp;
  Polyhedron *p1, *p2;

  /* Some basic checks: we only accept convex domains, with same constraint
   * and dimension numbers.
   */
  if (!cloog_domain_isconvex (d1) || !cloog_domain_isconvex (d2) ||
      (cloog_domain_nbconstraints (d1) != cloog_domain_nbconstraints (d2)) ||
      (cloog_domain_dim (d1) != cloog_domain_dim (d2)))
    return 0;

  p1 = d2p (d1);
  p2 = d2p (d2);
  nbc = cloog_domain_nbconstraints (d1);
  dim1 = cloog_domain_dim (d1);
  dim2 = cloog_domain_dim (d2);

  /* There should be only one difference between the two domains, it
   * has to be at the constant level and the difference must be of +1,
   * moreover, after the difference all domain coefficient has to be 0.
   * The matrix shape is:
   *
   * |===========|=====|<- 0 line
   * |===========|=====|
   * |===========|====?|<- different_constraint line (found here)
   * |===========|0000=|
   * |===========|0000=|<- pX->NbConstraints line
   *  ^         ^     ^
   *  |         |     |
   *  |         |     (pX->Dimension + 2) column
   *  |         scattdims column
   *  0 column
   */

  value_init_c (temp);
  for (i = 0; i < nbc; i++)
    {
      if (difference == 0)
	{			/* All elements except scalar must be equal. */
	  for (j = 0; j < dim1 + 1; j++)
	    if (value_ne (p1->Constraint[i][j],
			  p2->Constraint[i][j]))
	      {
		value_clear_c (temp);
		Polyhedron_Free (p1);
		Polyhedron_Free (p2);
		return 0;
	      }
	  /* The scalar may differ from +1 (now j=(p1->Dimension + 1)). */
	  if (value_ne (p1->Constraint[i][j],
			p2->Constraint[i][j]))
	    {
	      value_increment (temp, p2->Constraint[i][j]);
	      if (value_ne (p1->Constraint[i][j], temp))
		{
		  value_clear_c (temp);
		  Polyhedron_Free (p1);
		  Polyhedron_Free (p2);
		  return 0;
		}
	      else
		{
		  difference = 1;
		  different_constraint = i;
		}
	    }
	}
      else
	{			/* Scattering coefficients must be equal. */
	  for (j = 0; j < (scattdims + 1); j++)
	    if (value_ne (p1->Constraint[i][j],
			  p2->Constraint[i][j]))
	      {
		value_clear_c (temp);
		Polyhedron_Free (p1);
		Polyhedron_Free (p2);
		return 0;
	      }

	  /* Domain coefficients must be 0. */
	  for (; j < dim1 + 1; j++)
	    if (value_notzero_p (p1->Constraint[i][j])
		|| value_notzero_p (p2->Constraint[i][j]))
	      {
		value_clear_c (temp);
		Polyhedron_Free (p1);
		Polyhedron_Free (p2);
		return 0;
	      }

	  /* Scalar must be equal. */
	  if (value_ne (p1->Constraint[i][j],
			p2->Constraint[i][j]))
	    {
	      value_clear_c (temp);
	      Polyhedron_Free (p1);
	      Polyhedron_Free (p2);
	      return 0;
	    }
	}
    }
  value_clear_c (temp);

  /* If the domains are exactly the same, this is a block. */
  if (difference == 0)
    {
      Polyhedron_Free (p1);
      Polyhedron_Free (p2);
      return 1;
    }

  /* Now a basic check that the constraint with the difference is an
   * equality of a dimension with a constant.
   */
  for (i = 0; i <= different_constraint; i++)
    if (value_notzero_p (p1->Constraint[different_constraint][i]))
      {
	Polyhedron_Free (p1);
	Polyhedron_Free (p2);
	return 0;
      }

  if (value_notone_p (p1->Constraint[different_constraint][different_constraint + 1]))
    {
      Polyhedron_Free (p1);
      Polyhedron_Free (p2);
      return 0;
    }

  for (i = different_constraint + 2; i < dim1 + 1; i++)
    if (value_notzero_p (p1->Constraint[different_constraint][i]))
      {
	Polyhedron_Free (p1);
	Polyhedron_Free (p2);
	return 0;
      }

  /* For the moment, d1 and d2 are a block candidate. There remains to check
   * that there is no other domain that may put an integral point between
   * them. In our lazy test we ensure this property by verifying that the
   * constraint matrices have a very strict shape: let us consider that the
   * dimension with the difference is d. Then the first d dimensions are
   * defined in their depth order using equalities (thus the first column begins
   * with d zeroes, there is a d*d identity matrix and a zero-matrix for
   * the remaining simensions). If a domain can put integral points between the
   * domains of the block candidate, this means that the other entries on the
   * first d constraints are equal to those of d1 or d2. Thus we are looking for
   * such a constraint system, if it exists d1 and d2 is considered to not be
   * a block, it is a bock otherwise.
   *
   *  1. Only equalities (for the first different_constraint+1 lines).
   *  |  2. Must be the identity.
   *  |  |    3. Must be zero.
   *  |  |    |     4. Elements are equal, the last one is either date1 or 2.
   *  |  |    |     |
   *  | /-\ /---\ /---\
   * |0|100|00000|=====|<- 0 line
   * |0|010|00000|=====|
   * |0|001|00000|====?|<- different_constraint line
   * |*|***|*****|*****|
   * |*|***|*****|*****|<- pX->NbConstraints line
   *  ^   ^     ^     ^
   *  |   |     |     |
   *  |   |     |     (pX->Dimension + 2) column
   *  |   |     scattdims column
   *  |   different_constraint+1 column
   *  0 column
   */

  /* Step 1 and 2. This is only necessary to check one domain because
   * we checked that they are equal on this part before.
   */
  for (i = 0; i <= different_constraint; i++)
    {
      for (j = 0; j < i + 1; j++)
	if (value_notzero_p (p1->Constraint[i][j]))
	  {
	    Polyhedron_Free (p1);
	    Polyhedron_Free (p2);
	    return 0;
	  }

      if (value_notone_p (p1->Constraint[i][i + 1]))
	{
	  Polyhedron_Free (p1);
	  Polyhedron_Free (p2);
	  return 0;
	}

      for (j = i + 2; j <= different_constraint + 1; j++)
	if (value_notzero_p (p1->Constraint[i][j]))
	{
	  Polyhedron_Free (p1);
	  Polyhedron_Free (p2);
	  return 0;
	}
    }

  /* Step 3. */
  for (i = 0; i <= different_constraint; i++)
    for (j = different_constraint + 2; j <= scattdims; j++)
      if (value_notzero_p (p1->Constraint[i][j]))
	{
	  Polyhedron_Free (p1);
	  Polyhedron_Free (p2);
	  return 0;
	}

  value_init_c (date1);
  value_init_c (date2);
  value_init_c (date3);

  /* Now we have to check that the two different dates are unique. */
  value_assign (date1, p1->Constraint[different_constraint][dim1 + 1]);
  value_assign (date2, p2->Constraint[different_constraint][dim2 + 1]);

  /* Step 4. We check all domains except d1 and d2 and we look for at least
   * a difference with d1 or d2 on the first different_constraint+1 dimensions.
   */
  while (scattering != NULL)
    {
      if ((cloog_domain (scattering) != d1)
	  && (cloog_domain (scattering) != d2))
	{
	  CloogDomain *d3 = cloog_domain (scattering);
	  Polyhedron *p3 = d2p (d3);
	  int dim3 = cloog_domain_dim (d3);

	  value_assign (date3,
			p3->Constraint[different_constraint][dim3 + 1]);
	  difference = 0;

	  if (value_ne (date3, date2) && value_ne (date3, date1))
	    difference = 1;

	  for (i = 0; (i < different_constraint) && (!difference); i++)
	    for (j = 0; (j < dim3 + 2) && !difference; j++)
	      if (value_ne
		  (p1->Constraint[i][j],
		   p3->Constraint[i][j]))
		difference = 1;

	  for (j = 0; (j < dim3 + 1) && !difference; j++)
	    if (value_ne
		(p1->Constraint[different_constraint][j],
		 p3->Constraint[different_constraint][j]))
	      difference = 1;

	  Polyhedron_Free (p3);
	  if (!difference)
	    {
	      value_clear_c (date1);
	      value_clear_c (date2);
	      value_clear_c (date3);
	      Polyhedron_Free (p1);
	      Polyhedron_Free (p2);
	      return 0;
	    }
	}

      scattering = cloog_next_domain (scattering);
    }

  Polyhedron_Free (p1);
  Polyhedron_Free (p2);
  value_clear_c (date1);
  value_clear_c (date2);
  value_clear_c (date3);
  return 1;
}


/**
 * cloog_domain_lazy_disjoint function:
 * This function returns 1 if the domains given as input are disjoint, 0 if it
 * is unable to decide. This function finds the unknown with fixed values in
 * both domains (on a given constraint, their column entry is not zero and
 * only the constant coefficient can be different from zero) and verify that
 * their values are the same. If not, the domains are obviously disjoint and
 * it returns 1, if there is not such case it returns 0.  This is a very fast
 * way to verify this property. It has been shown (with the CLooG benchmarks)
 * that operations on disjoint domains are 36% of all the polyhedral
 * computations. For 94% of the actually identical domains, this
 * function answer that they are disjoint and allow to give immediately the
 * trivial solution instead of calling the heavy general functions of PolyLib.
 * - August 22th 2003: first version.
 * - June   21rd 2005: Adaptation for GMP (based on S. Verdoolaege's version of
 *                     CLooG 0.12.1).
 */
int
cloog_domain_lazy_disjoint (CloogDomain * d1, CloogDomain * d2)
{
  int i1, j1, i2, j2, scat_dim, nbc1, nbc2;
  int dim1, dim2;
  Value scat_val;
  Polyhedron *p1, *p2;

  if (!cloog_domain_isconvex (d1) || !cloog_domain_isconvex (d2))
    return 0;

  p1 = d2p (d1);
  p2 = d2p (d2);
  nbc1 = cloog_domain_nbconstraints (d1);
  nbc2 = cloog_domain_nbconstraints (d2);
  dim1 = cloog_domain_dim (d1);
  dim2 = cloog_domain_dim (d2);
  value_init_c (scat_val);

  for (i1 = 0; i1 < nbc1; i1++)
    {
      if (value_notzero_p (p1->Constraint[i1][0]))
	continue;

      scat_dim = 1;
      while (value_zero_p (p1->Constraint[i1][scat_dim]) &&
	     (scat_dim < dim1))
	scat_dim++;

      if (value_notone_p (p1->Constraint[i1][scat_dim]))
	continue;
      else
	{
	  for (j1 = scat_dim + 1; j1 <= dim1; j1++)
	    if (value_notzero_p (p1->Constraint[i1][j1]))
	      break;

	  if (j1 != dim1 + 1)
	    continue;

	  value_assign (scat_val,
			p1->Constraint[i1][dim1 + 1]);

	  for (i2 = 0; i2 < nbc2; i2++)
	    {
	      for (j2 = 0; j2 < scat_dim; j2++)
		if (value_notzero_p (p2->Constraint[i2][j2]))
		  break;

	      if ((j2 != scat_dim)
		  ||
		  value_notone_p (p2->Constraint[i2][scat_dim]))
		continue;

	      for (j2 = scat_dim + 1; j2 < dim2; j2++)
		if (value_notzero_p (p2->Constraint[i2][j2]))
		  break;

	      if (j2 != dim2)
		continue;

	      if (value_ne
		  (p2->Constraint[i2][dim2 + 1], scat_val))
		{
		  value_clear_c (scat_val);
		  return 1;
		}
	    }
	}
    }

  value_clear_c (scat_val);
  Polyhedron_Free (p1);
  Polyhedron_Free (p2);
  return 0;
}


/**
 * cloog_domain_list_lazy_same function:
 * This function returns 1 if two domains in the list are the same, 0 if it
 * is unable to decide.
 * - February 9th 2004: first version.
 */
int
cloog_domain_list_lazy_same (CloogDomainList * list)
{				/*int i=1, j=1 ; */
  CloogDomainList *current, *next;

  current = list;
  while (current != NULL)
    {
      next = cloog_next_domain (current);
      /*j=i+1; */
      while (next != NULL)
	{
	  if (cloog_domain_lazy_equal (cloog_domain (current),
				       cloog_domain (next)))
	    {			/*printf("Same domains: %d and %d\n",i,j) ; */
	      return 1;
	    }
	  /*j++ ; */
	  next = cloog_next_domain (next);
	}
      /*i++ ; */
      current = cloog_next_domain (current);
    }

  return 0;
}

/**
 * cloog_domain_cut_first function:
 * this function returns a CloogDomain structure with everything except the
 * first part of the polyhedra union of the input domain as domain. After a call
 * to this function, there remains in the CloogDomain structure provided as
 * input only the first part of the original polyhedra union.
 * - April 20th 2005: first version, extracted from different part of loop.c.
 */
CloogDomain *
cloog_domain_cut_first (CloogDomain * domain)
{
  CloogDomain *rest;

  if (domain && cloog_domain_polyhedron (domain))
    {
      if (!cloog_upol_next (cloog_domain_upol (domain)))
	return NULL;

      rest = cloog_new_domain (cloog_upol_next (cloog_domain_upol (domain)));
      cloog_upol_set_next (cloog_domain_upol (domain), NULL);
    }
  else
    rest = NULL;

  return print_result ("cloog_domain_cut_first", cloog_check_domain (rest));
}


/**
 * cloog_domain_lazy_isscalar function:
 * this function returns 1 if the dimension 'dimension' in the domain 'domain'
 * is scalar, this means that the only constraint on this dimension must have
 * the shape "x.dimension + scalar = 0" with x an integral variable. This
 * function is lazy since we only accept x=1 (further calculations are easier
 * in this way).
 * - June 14th 2005: first version.
 * - June 21rd 2005: Adaptation for GMP.
 */
int
cloog_domain_lazy_isscalar (CloogDomain * domain, int dimension)
{
  int i, j;
  Polyhedron *polyhedron = d2p (domain);
  int nbc = cloog_domain_nbconstraints (domain);
  int dim = cloog_domain_dim (domain);

  /* For each constraint... */
  for (i = 0; i < nbc; i++)
    {				/* ...if it is concerned by the potentially scalar dimension... */
      if (value_notzero_p
	  (polyhedron->Constraint[i][dimension + 1]))
	{			/* ...check that the constraint has the shape "dimension + scalar = 0". */
	  for (j = 0; j <= dimension; j++)
	    if (value_notzero_p (polyhedron->Constraint[i][j]))
	      {
		Polyhedron_Free (polyhedron);
		return 0;
	      }

	  if (value_notone_p
	      (polyhedron->Constraint[i][dimension + 1]))
	    {
	      Polyhedron_Free (polyhedron);
	      return 0;
	    }

	  for (j = dimension + 2; j < dim + 1; j++)
	    if (value_notzero_p (polyhedron->Constraint[i][j]))
	      {
		Polyhedron_Free (polyhedron);
		return 0;
	      }
	}
    }

  Polyhedron_Free (polyhedron);
  return 1;
}


/**
 * cloog_domain_scalar function:
 * when we call this function, we know that "dimension" is a scalar dimension,
 * this function finds the scalar value in "domain" and returns it in "value".
 * - June 30th 2005: first version.
 */
void
cloog_domain_scalar (CloogDomain * domain, int dimension, Value * value)
{
  int i;
  Polyhedron *polyhedron = d2p (domain);
  int nbc = cloog_domain_nbconstraints (domain);
  int dim = cloog_domain_dim (domain);

  /* For each constraint... */
  for (i = 0; i < nbc; i++)
    {				/* ...if it is the equality defining the scalar dimension... */
      if (value_notzero_p
	  (polyhedron->Constraint[i][dimension + 1])
	  && value_zero_p (polyhedron->Constraint[i][0]))
	{			/* ...Then send the scalar value. */
	  value_assign (*value, polyhedron->Constraint[i][dim + 1]);
	  value_oppose (*value, *value);
	  Polyhedron_Free (polyhedron);
	  return;
	}
    }

  /* We should have found a scalar value: if not, there is an error. */
  fprintf (stderr, "[CLooG]ERROR: dimension %d is not scalar as expected.\n",
	   dimension);
  Polyhedron_Free (polyhedron);
  exit (0);
}


/**
 * cloog_domain_erase_dimension function:
 * this function returns a CloogDomain structure builds from 'domain' where
 * we removed the dimension 'dimension' and every constraint considering this
 * dimension. This is not a projection ! Every data concerning the
 * considered dimension is simply erased.
 * - June 14th 2005: first version.
 * - June 21rd 2005: Adaptation for GMP.
 */
CloogDomain *
cloog_domain_erase_dimension (CloogDomain * domain, int dimension)
{
  int i, j, mi, nb_dim, nbc;
  CloogMatrix *matrix;
  CloogDomain *erased;
  Polyhedron *polyhedron;

  polyhedron = d2p (domain);
  nb_dim = cloog_domain_dim (domain);
  nbc = cloog_domain_nbconstraints (domain);

  /* The matrix is one column less and at least one constraint less. */
  matrix = cloog_matrix_alloc (nbc - 1, nb_dim + 1);

  /* mi is the constraint counter for the matrix. */
  mi = 0;
  for (i = 0; i < nbc; i++)
    if (value_zero_p (polyhedron->Constraint[i][dimension + 1]))
      {
	for (j = 0; j <= dimension; j++)
	  value_assign (matrix->p[mi][j],
			polyhedron->Constraint[i][j]);

	for (j = dimension + 2; j < nb_dim + 2; j++)
	  value_assign (matrix->p[mi][j - 1],
			polyhedron->Constraint[i][j]);

	mi++;
      }

  erased = cloog_domain_matrix2domain (matrix);
  cloog_matrix_free (matrix);

  Polyhedron_Free (polyhedron);
  return print_result ("cloog_domain_erase_dimension", cloog_check_domain (erased));
}

/* Number of polyhedra inside the union of disjoint polyhedra.  */

unsigned
cloog_domain_nb_polyhedra (CloogDomain * domain)
{
  unsigned res = 0;
  ppl_polyhedra_union *upol = cloog_domain_upol (domain);

  while (upol)
    {
      res++;
      upol = cloog_upol_next (upol);
    }

  return res;
}

void
cloog_domain_print_polyhedra (FILE * foo, CloogDomain * domain)
{
  ppl_polyhedra_union *upol = cloog_domain_upol (domain);

  while (upol != NULL)
    {
      CloogMatrix *matrix = cloog_upol_domain2matrix (upol);
      cloog_matrix_print (foo, matrix);
      cloog_matrix_free (matrix);
      upol = cloog_upol_next (upol);
    }
}

void
debug_cloog_domain (CloogDomain *domain)
{
  cloog_domain_print_polyhedra (stderr, domain);
}

void
debug_cloog_matrix (CloogMatrix *m)
{
  cloog_matrix_print (stderr, m);
}
