
   /**-------------------------------------------------------------------**
    **                               CLooG                               **
    **-------------------------------------------------------------------**
    **                             checking.c                            **
    **-------------------------------------------------------------------**
    **                   First version: July 22th 2008                   **
    **-------------------------------------------------------------------**/


/******************************************************************************
 *               CLooG : the Chunky Loop Generator (experimental)             *
 ******************************************************************************
 *                                                                            *
 * Copyright (C) 2001-2008 Cedric Bastoul                                     *
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

static CloogDomain *
cloog_domain_addconstraints_1 (domain_source, domain_target)
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

CloogDomain *
cloog_domain_extend_1 (CloogDomain * domain, int dim, int nb_par)
{
  int row, column, nb_rows, nb_columns, difference;
  CloogDomain *extended_domain;
  CloogMatrix *matrix;

  nb_rows = 1 + cloog_domain_dim (domain);
  nb_columns = dim + nb_par + 1;
  difference = nb_columns - nb_rows;

  if (difference == 0)
    return print_result ("cloog_domain_extend_1", cloog_domain_copy (domain));

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

  return print_result ("cloog_domain_extend_1", cloog_check_domain (extended_domain));
}

static CloogDomain *
cloog_domain_convex_1 (CloogDomain * domain)
{
  Polyhedron *p = d2p (domain);
  CloogDomain *res =
    cloog_check_domain (cloog_domain_alloc
			(DomainConvex (p, MAX_RAYS)));
  Polyhedron_Free (p);
  return print_result ("cloog_domain_convex_1", res);
}

static CloogDomain *
cloog_domain_simplify_1 (CloogDomain * dom1, CloogDomain * dom2)
{
  CloogMatrix *M, *M2;
  CloogDomain *dom;
  Polyhedron *p1 = d2p (dom1);
  Polyhedron *p2 = d2p (dom2);
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
	  Vector_Copy (p1->Constraint[i], M->p[rank], cols);
	  if (Gauss (M, rank + 1, cols - 1) > rank)
	    {
	      Vector_Copy (p1->Constraint[i], M2->p[row++], cols);
	      rank++;
	    }
	}
      if (row < cloog_domain_nbeq (dom1))
	{
	  Vector_Copy (p1->Constraint[cloog_domain_nbeq (dom1)],
		       M2->p[row], (nbc - cloog_domain_nbeq (dom1)) * cols);
	  p1 = Constraints2Polyhedron (M2, MAX_RAYS);
	}
      cloog_matrix_free (M2);
      cloog_matrix_free (M);
    }

  dom = cloog_domain_alloc (DomainSimplify (p1, p2, MAX_RAYS));
  Polyhedron_Free (p1);
  Polyhedron_Free (p2);
  return print_result ("cloog_domain_simplify", cloog_check_domain (dom));
}

