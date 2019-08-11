/***************************************************************************
 *            hybrid_set.cc
 *
 *  Copyright 2008  Pieter Collins
 * 
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */
 
#include "hybrid_set.h"
#include "hybrid_automaton_interface.h"

namespace Ariadne {


/* HybridBox *************************************************************************************/


HybridBoxes::
HybridBoxes(
		const HybridSpace& space,
		const Box& bbox)
{
    for(HybridSpace::const_iterator loc_iter = space.begin(); loc_iter != space.end(); ++loc_iter) {
    	ARIADNE_ASSERT_MSG(loc_iter->second == bbox.dimension(),
    			"The provided box dimension is not the same as that of location " << loc_iter->first.name());
        this->insert(make_pair(loc_iter->first,bbox));
    }
}

HybridBoxes::
HybridBoxes(const HybridSpace& space)
{
    for(HybridSpace::const_iterator loc_iter = space.begin(); loc_iter != space.end(); ++loc_iter)
        this->insert(make_pair(loc_iter->first,Box::empty_box(loc_iter->second)));
}


tribool
HybridBoxes::
overlaps(const LocalisedBox& hbx) const
{
	locations_const_iterator loc_iter=this->find(hbx.first);
	if (loc_iter!=this->locations_end()) return loc_iter->second.overlaps(hbx.second);
	else return false;
}


tribool
HybridBoxes::
disjoint(const LocalisedBox& hbx) const
{
	locations_const_iterator loc_iter=this->find(hbx.first);
	if (loc_iter!=this->locations_end()) return loc_iter->second.disjoint(hbx.second);
	else return true;
}


tribool
HybridBoxes::
covers(const LocalisedBox& hbx) const
{
	locations_const_iterator loc_iter=this->find(hbx.first);
	if (loc_iter!=this->locations_end()) return loc_iter->second.covers(hbx.second);
	else return false;
}


tribool
HybridBoxes::
inside(const HybridBoxes& hbx) const
{
	ARIADNE_ASSERT_MSG(this->space() == hbx.space(), "The two HybridBox have different spaces.");

	tribool result = true;
    for(locations_const_iterator loc_iter=hbx.begin(); loc_iter!=hbx.end(); ++loc_iter) {
    	locations_const_iterator this_iter = this->find(loc_iter->first);

    	tribool local_result = this_iter->second.inside(loc_iter->second);
    	if (!possibly(local_result))
    		return false;
    	else if (indeterminate(local_result))
    		result = indeterminate;
	}
	return result;
}


tribool
HybridBoxes::
superset(const HybridBoxes& hbx) const
{
	tribool result = true;

	ARIADNE_ASSERT_MSG(this->space() == hbx.space(), "The boxes must have the same space.");
	for (locations_const_iterator loc_it = this->begin(); loc_it != this->end(); ++loc_it) {
		tribool local_result = loc_it->second.superset(hbx.find(loc_it->first)->second);
		if (definitely(!local_result))
			return false;
		else if (indeterminate(local_result))
			result = indeterminate;
	}
	return result;
}




Box
HybridBoxes::
project(const std::vector<uint>& dimensions) const
{
	ARIADNE_ASSERT_MSG(dimensions.size()>0, "Provide at least one dimension to project to.");

	Box result = Box::empty_box(dimensions.size());

	for (locations_const_iterator loc_it = this->begin(); loc_it != this->end(); ++loc_it)
		result = hull(result,loc_it->second.project(dimensions));

	return result;
}


/* HybridGrid *********************************************************************************/


std::map<DiscreteLocation,Vector<Float> >
HybridGrid::
lengths() const
{
	std::map<DiscreteLocation,Vector<Float> > result;

	for (HybridGrid::const_iterator it = this->begin(); it != this->end(); ++it)
		result.insert(make_pair(it->first,it->second.lengths()));

	return result;
}


/* HybridConstraintSet ************************************************************************/


HybridConstraintSet::
HybridConstraintSet(const HybridVectorFunction& func,const HybridBoxes& codomain)
{
	ARIADNE_ASSERT_MSG(codomain.size() == func.size(), "The codomain and functions have different sizes.");

	for (HybridBoxes::const_iterator cod_it = codomain.begin(); cod_it != codomain.end(); ++cod_it) {
		HybridVectorFunction::const_iterator func_it = func.find(cod_it->first);
		ARIADNE_ASSERT_MSG(func_it != func.end(), "No function for location " << cod_it->first.name() << " is present.");
		this->insert(make_pair(cod_it->first,ConstraintSet(func_it->second,cod_it->second)));
	}
}


HybridConstraintSet::
HybridConstraintSet(const HybridSpace& hspace, ConstraintSet constraint)
{
	for (HybridSpace::const_iterator hs_it = hspace.begin(); hs_it != hspace.end(); ++hs_it) {
		ARIADNE_ASSERT_MSG(hs_it->second == constraint.function().argument_size(),
				"The continuous space would not match the constraint function.");
		this->insert(make_pair(hs_it->first,constraint));
	}
}


HybridConstraintSet::
HybridConstraintSet(const HybridSpace& hspace)
{
	for (HybridSpace::const_iterator hs_it = hspace.begin(); hs_it != hspace.end(); ++hs_it) {
		this->insert(make_pair(hs_it->first,ConstraintSet(hs_it->second)));
	}
}


/* HybridBoundedConstraintSet **************************************************************************/


HybridBoundedConstraintSet::
HybridBoundedConstraintSet(
		const HybridBoxes& domain,
		const HybridVectorFunction& func,
		const HybridBoxes& codomain)
{
	ARIADNE_ASSERT_MSG(domain.size() == func.size(), "The domain and functions have different sizes.");
	ARIADNE_ASSERT_MSG(codomain.size() == func.size(), "The codomain and functions have different sizes.");

	for (HybridBoxes::const_iterator cod_it = codomain.begin(); cod_it != codomain.end(); ++cod_it) {
		HybridBoxes::const_iterator dom_it = domain.find(cod_it->first);
		ARIADNE_ASSERT_MSG(dom_it != domain.end(), "No domain for location " << cod_it->first.name() << " is present.");
		HybridVectorFunction::const_iterator func_it = func.find(cod_it->first);
		ARIADNE_ASSERT_MSG(func_it != func.end(), "No function for location " << cod_it->first.name() << " is present.");

		this->insert(make_pair(cod_it->first,BoundedConstraintSet(dom_it->second,func_it->second,cod_it->second)));
	}
}


HybridBoundedConstraintSet::
HybridBoundedConstraintSet(
		const HybridSpace& hspace,
		BoundedConstraintSet constraint)
{
	for (HybridSpace::const_iterator hs_it = hspace.begin(); hs_it != hspace.end(); ++hs_it) {
		ARIADNE_ASSERT_MSG(hs_it->second == constraint.function().argument_size(),
				"The continuous space would not match the constraint function.");
		this->insert(make_pair(hs_it->first,constraint));
	}
}


HybridBoundedConstraintSet::
HybridBoundedConstraintSet(const HybridSpace& hspace, const Box& domain_bx)
{
	for (HybridSpace::const_iterator hs_it = hspace.begin(); hs_it != hspace.end(); ++hs_it) {
		ARIADNE_ASSERT_MSG(hs_it->second == domain_bx.dimension(),
				"The box dimension does not match location " << hs_it->first.name() << " dimension.");
		this->insert(make_pair(hs_it->first,BoundedConstraintSet(domain_bx)));
	}
}

HybridBoundedConstraintSet::
HybridBoundedConstraintSet(const HybridSpace& hspace)
{
	for (HybridSpace::const_iterator hs_it = hspace.begin(); hs_it != hspace.end(); ++hs_it) {
		this->insert(make_pair(hs_it->first,BoundedConstraintSet(Box::empty_box(hs_it->second))));
	}
}


HybridBoundedConstraintSet::
HybridBoundedConstraintSet(const HybridBoxes& domain_boxes)
{
	for (HybridBoxes::const_iterator hbx_it = domain_boxes.begin(); hbx_it != domain_boxes.end(); ++hbx_it) {
		this->insert(make_pair(hbx_it->first,BoundedConstraintSet(hbx_it->second)));
	}
}


/* Free functions ********************************************************************************/


HybridBoxes
hull(const HybridBoxes& box1, const HybridBoxes& box2)
{
	ARIADNE_ASSERT_MSG(box1.size() == box2.size(),"The two hybrid boxes must have the same number of locations ("
												  << box1.size() << " vs " << box2.size() << ").");
	HybridBoxes result;
	for (HybridBoxes::const_iterator box1_it = box1.begin(); box1_it != box1.end(); ++box1_it) {
		HybridBoxes::const_iterator box2_it = box2.find(box1_it->first);
		ARIADNE_ASSERT_MSG(box2_it != box2.end(),"The location " << box1_it->first.name() << " is not present in both hybrid boxes.");
		result.insert(std::pair<DiscreteLocation,Box>(box1_it->first,hull(box1_it->second,box2_it->second)));
	}

	return result;
}


HybridBoxes
widen(const HybridBoxes& box, const HybridFloatVector& epsilon)
{
	HybridBoxes result;

	for (HybridBoxes::const_iterator loc_it = box.begin(); loc_it != box.end(); ++loc_it) {
		HybridFloatVector::const_iterator epsilon_it = epsilon.find(loc_it->first);
		ARIADNE_ASSERT_MSG(epsilon_it != epsilon.end(),"The location " << loc_it->first.name() << " is not present in the epsilon map.");
		Box widened_box = loc_it->second;
		widened_box.widen(epsilon_it->second);
		result.insert(std::pair<DiscreteLocation,Box>(loc_it->first,widened_box));
	}

	return result;
}


HybridBoxes
widen(const HybridBoxes& box)
{
	HybridBoxes result;
	for (HybridBoxes::const_iterator hbx_it = box.begin(); hbx_it != box.end(); hbx_it++) {	
                DiscreteLocation loc = hbx_it->first;
                Box bx = hbx_it->second;
		bx.widen();
		result.insert(make_pair(loc,bx));
	}

	return result;
}


bool subset(const HybridDenotableSet& theSet1, const HybridDenotableSet& theSet2)
{
	ARIADNE_ASSERT_MSG(theSet1.grid() == theSet2.grid(), "Cannot perform subset on two HybridDenotableSets with different grids.")
	for (HybridDenotableSet::locations_const_iterator set1_it = theSet1.locations_begin();
			set1_it != theSet1.locations_end(); ++set1_it) {
		if (!subset(set1_it->second,theSet2.find(set1_it->first)->second))
			return false;
	}

	return true;
}

tribool disjoint(const HybridConstraintSet& cons_set, const HybridDenotableSet& grid_set)
{
	ARIADNE_ASSERT_MSG(cons_set.space() == grid_set.space(), "The denotable set and constraint set have mismatched spaces.");

    tribool result = true;

    for (HybridDenotableSet::locations_const_iterator gts_it = grid_set.locations_begin(); gts_it != grid_set.locations_end(); ++gts_it) {
    	HybridConstraintSet::locations_const_iterator cs_it = cons_set.find(gts_it->first);
    	if (cs_it != cons_set.locations_end()) {
    		tribool disjoint_gts = disjoint(cs_it->second,gts_it->second);
    		if (!disjoint_gts)
    			return false;
    		else
    			result = result && disjoint_gts;
    	}
    }
    return result;
}


tribool overlaps(const HybridConstraintSet& cons_set, const HybridDenotableSet& grid_set)
{
    return !disjoint(cons_set,grid_set);
}


tribool covers(const HybridConstraintSet& cons_set, const HybridDenotableSet& grid_set)
{
	ARIADNE_ASSERT_MSG(cons_set.space() == grid_set.space(), "The denotable set and constraint set have mismatched spaces.");

    tribool result = true;

    for (HybridDenotableSet::locations_const_iterator gts_it = grid_set.locations_begin(); gts_it != grid_set.locations_end(); ++gts_it) {
    	HybridConstraintSet::locations_const_iterator cs_it = cons_set.find(gts_it->first);
    	if (cs_it != cons_set.locations_end()) {
    		tribool covered_gts = covers(cs_it->second,gts_it->second);
    		if (!covered_gts)
    			return false;
    		else
    			result = result && covered_gts;
    	}
    }
    return result;

}


HybridDenotableSet outer_intersection(const HybridDenotableSet& hds_set, const HybridConstraintSet& cons_set)
{
	HybridDenotableSet result;

    for (HybridDenotableSet::locations_const_iterator hds_it = hds_set.locations_begin(); hds_it != hds_set.locations_end(); ++hds_it) {
    	HybridConstraintSet::locations_const_iterator cs_it = cons_set.find(hds_it->first);
    	if (cs_it != cons_set.locations_end())
    		result.insert(make_pair(hds_it->first,outer_intersection(hds_it->second,cs_it->second)));
    }

	return result;
}


HybridDenotableSet inner_intersection(const HybridDenotableSet& hds_set, const HybridConstraintSet& cons_set)
{
	HybridDenotableSet result;

    for (HybridDenotableSet::locations_const_iterator hds_it = hds_set.locations_begin(); hds_it != hds_set.locations_end(); ++hds_it) {
    	HybridConstraintSet::locations_const_iterator cs_it = cons_set.find(hds_it->first);
    	if (cs_it != cons_set.locations_end())
    		result.insert(make_pair(hds_it->first,inner_intersection(hds_it->second,cs_it->second)));
    }

	return result;
}


HybridDenotableSet outer_difference(const HybridDenotableSet& hd_set, const HybridConstraintSet& cons_set)
{
	HybridDenotableSet result;

	for (HybridDenotableSet::locations_const_iterator hds_it = hd_set.locations_begin();
			hds_it != hd_set.locations_end(); ++hds_it) {
		const DiscreteLocation& loc = hds_it->first;
		const DenotableSetType& ds = hds_it->second;

		if (ds.empty()) {
			result.insert(*hds_it);
		} else {
			DenotableSetType local_result = ds;
			HybridConstraintSet::const_iterator cons_it = cons_set.find(loc);
		    if (cons_it != cons_set.end()) {
		    	DenotableSetType inner_intersection = Ariadne::inner_intersection(ds,cons_it->second);
		    	local_result.remove(inner_intersection);
		    }
			result.insert(make_pair(loc,local_result));
		}
	}
	return result;
}


HybridDenotableSet inner_difference(const HybridDenotableSet& hd_set, const HybridConstraintSet& cons_set)
{
	HybridDenotableSet result;

	for (HybridDenotableSet::locations_const_iterator hds_it = hd_set.locations_begin();
			hds_it != hd_set.locations_end(); ++hds_it) {
		const DiscreteLocation& loc = hds_it->first;
		const DenotableSetType& ds = hds_it->second;

		if (ds.empty()) {
			result.insert(*hds_it);
		} else {
			DenotableSetType local_result = ds;
			HybridConstraintSet::const_iterator cons_it = cons_set.find(loc);
		    if (cons_it != cons_set.end()) {
		    	DenotableSetType outer_intersection = Ariadne::outer_intersection(ds,cons_it->second);
		    	local_result.remove(outer_intersection);
		    }
			result.insert(make_pair(loc,local_result));
		}
	}
	return result;
}


HybridBoxes eps_codomain(
		const HybridDenotableSet& grid_set,
		const HybridFloatVector& eps,
		const HybridVectorFunction& func)
{
	HybridBoxes result;

	for (HybridDenotableSet::locations_const_iterator loc_it = grid_set.locations_begin(); loc_it != grid_set.locations_end(); ++loc_it) {
		const DiscreteLocation& loc = loc_it->first;
		HybridFloatVector::const_iterator eps_it = eps.find(loc);
		HybridVectorFunction::const_iterator func_it = func.find(loc);
		ARIADNE_ASSERT(eps_it != eps.end());
		if (func_it != func.end())
			result.insert(std::pair<DiscreteLocation,Box>(loc,eps_codomain(loc_it->second,eps_it->second,func_it->second)));
	}

	return result;
}


DenotableSetType flatten_and_project_down(const HybridDenotableSet& grid_set, const Vector<uint>& indices)
{
	const HybridGrid hgrid = grid_set.grid();

	HybridGrid::const_iterator hg_it = hgrid.begin();
	const Grid& grid = hg_it->second;
	++hg_it;
	for (; hg_it != hgrid.end(); ++hg_it) {
		ARIADNE_ASSERT_MSG(hg_it->second == grid,"The grid must be the same for all locations.");
	}

	Grid projected_grid = project_down(grid,indices);
	DenotableSetType result(projected_grid);

	for (HybridDenotableSet::locations_const_iterator gts_it = grid_set.locations_begin();
			gts_it != grid_set.locations_end(); ++gts_it)
		result.adjoin(project_down(gts_it->second,indices));

	return result;
}


HybridFloatVector
getHybridMidpointAbsoluteDerivatives(
		const HybridAutomatonInterface& sys,
		const HybridBoxes& bounding_domain)
{
	HybridFloatVector result;

	const HybridSpace hspace = sys.state_space();
	for (HybridSpace::const_iterator hs_it = hspace.begin(); hs_it != hspace.end(); hs_it++) {

		const DiscreteLocation& loc = hs_it->first;
		const uint dim = hs_it->second;
		HybridBoxes::const_iterator domain_box_it = bounding_domain.find(loc);

		ARIADNE_ASSERT_MSG(domain_box_it != bounding_domain.end(),
				"The system state space and the domain space do not match.");

		const Box& domain_box = domain_box_it->second;

		Vector<Float> local_result(dim,0);

		if (!domain_box.empty()) {
			Vector<Interval> der_bbx = sys.dynamic_function(loc)(domain_box);
			for (uint i=0;i<dim;i++)
				local_result[i] = abs(der_bbx[i]).midpoint();
		}

		result.insert(make_pair(loc,local_result));
	}

	return result;
}


} // namespace Ariadne
