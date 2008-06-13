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
