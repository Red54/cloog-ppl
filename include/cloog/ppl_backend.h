/* Copyright (C) 2008 Sebastian Pop                                        

   This is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This software is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with software; if not, write to the Free Software Foundation,
   Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
*/

# include <ppl_c.h>


extern void cloog_initialize (void);

static inline void
cloog_finalize (void)
{
  ppl_finalize ();
}

typedef struct ppl_polyhedra_union {
  Polyhedron *_polyhedron;
  struct ppl_polyhedra_union *_next;
} ppl_polyhedra_union;

typedef struct cloogdomain
{
  struct ppl_polyhedra_union *_polyhedron;
  int _references;
} CloogDomain;

extern void debug_cloog_domain (CloogDomain *);

static inline ppl_polyhedra_union *
cloog_domain_upol (CloogDomain * domain)
{
  return domain->_polyhedron;
}

static inline Polyhedron *
cloog_domain_polyhedron (CloogDomain * domain)
{
  return domain->_polyhedron->_polyhedron;
}

static inline Polyhedron *
cloog_upol_polyhedron (ppl_polyhedra_union * ppl)
{
  return ppl->_polyhedron;
}

static inline void
cloog_upol_set_polyhedron (ppl_polyhedra_union * ppl, Polyhedron *p)
{
  ppl->_polyhedron = p;
}

static inline ppl_polyhedra_union *
cloog_new_upol (Polyhedron *p)
{
  ppl_polyhedra_union *ppl = (ppl_polyhedra_union *) malloc (sizeof (ppl_polyhedra_union));
  ppl->_polyhedron = p;
  ppl->_next = NULL;
  return ppl;
}

static inline void
cloog_upol_set_next (ppl_polyhedra_union * p, ppl_polyhedra_union * n)
{
  p->_next = n;
}

static inline ppl_polyhedra_union *
cloog_upol_next (ppl_polyhedra_union * p)
{
  return p->_next;
}


