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

static inline Polyhedron *
cloog_domain_polyhedron (CloogDomain * domain)
{
  return domain->_polyhedron->_polyhedron;
}
