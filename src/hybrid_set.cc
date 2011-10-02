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

namespace Ariadne {

std::map<DiscreteLocation,Vector<Float> >
HybridGrid::
lengths() const
{
	std::map<DiscreteLocation,Vector<Float> > result;

	for (HybridGrid::const_iterator it = this->begin(); it != this->end(); ++it)
		result.insert(std::pair<DiscreteLocation,Vector<Float> >(it->first,it->second.lengths()));

	return result;
}

HybridConstraintSet::
HybridConstraintSet(const HybridVectorFunction& func,const HybridBoxes& codomain)
{
	ARIADNE_ASSERT_MSG(codomain.size() == func.size(), "The codomain and functions have different sizes.");

	for (HybridBoxes::const_iterator cod_it = codomain.begin(); cod_it != codomain.end(); ++cod_it) {
		HybridVectorFunction::const_iterator func_it = func.find(cod_it->first);
		ARIADNE_ASSERT_MSG(func_it != func.end(), "No function for location " << cod_it->first.name() << " is present.");
		this->insert(std::pair<DiscreteLocation,ConstraintSet>(
				cod_it->first,ConstraintSet(func_it->second,cod_it->second)));
	}
}

HybridConstraintSet::
HybridConstraintSet(const HybridSpace& hspace, ConstraintSet constraint)
{
	for (HybridSpace::const_iterator hs_it = hspace.begin(); hs_it != hspace.end(); ++hs_it) {
		ARIADNE_ASSERT_MSG(hs_it->second == constraint.function().argument_size(),
				"The continuous space would not match the constraint function.");
		this->insert(std::pair<DiscreteLocation,ConstraintSet>(hs_it->first,constraint));
	}
}


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


bool
superset(const HybridBoxes& box1, const HybridBoxes& box2)
{
	for (HybridBoxes::const_iterator box1_it = box1.begin(); box1_it != box1.end(); ++box1_it) {
		HybridBoxes::const_iterator box2_it = box2.find(box1_it->first);
		ARIADNE_ASSERT_MSG(box2_it != box2.end(),"The location " << box1_it->first.name() << " is not present in both hybrid boxes.");
		if (!box1_it->second.superset(box2_it->second))
			return false;
	}

	return true;
}


bool
covers(const HybridBoxes& hboxes, const HybridBox& hbox)
{
	HybridBoxes::const_iterator hboxes_in_location = hboxes.find(hbox.first);
	ARIADNE_ASSERT_MSG(hboxes_in_location != hboxes.end(),
			"The location " << hbox.first << " is not present in the HybridBoxes.");

	return hboxes_in_location->second.covers(hbox.second);
}


HybridBoxes
shrink_in(const HybridBoxes& box, const HybridFloatVector& epsilon)
{
	HybridBoxes result;

	for (HybridBoxes::const_iterator loc_it = box.begin(); loc_it != box.end(); ++loc_it) {
		HybridFloatVector::const_iterator epsilon_it = epsilon.find(loc_it->first);
		ARIADNE_ASSERT_MSG(epsilon_it != epsilon.end(),"The location " << loc_it->first.name() << " is not present in the epsilon map.");
		result.insert(std::pair<DiscreteLocation,Box>(loc_it->first,loc_it->second.shrink_in(epsilon_it->second)));
	}

	return result;
}


HybridBoxes
shrink_out(const HybridBoxes& box, const HybridFloatVector& epsilon)
{
	HybridBoxes result;

	for (HybridBoxes::const_iterator loc_it = box.begin(); loc_it != box.end(); ++loc_it) {
		HybridFloatVector::const_iterator epsilon_it = epsilon.find(loc_it->first);
		ARIADNE_ASSERT_MSG(epsilon_it != epsilon.end(),"The location " << loc_it->first.name() << " is not present in the epsilon map.");
		result.insert(std::pair<DiscreteLocation,Box>(loc_it->first,loc_it->second.shrink_out(epsilon_it->second)));
	}

	return result;
}


HybridBoxes
widen(const HybridBoxes& box)
{
	HybridBoxes result;
	for (HybridBoxes::const_iterator loc_it = box.begin(); loc_it != box.end(); loc_it++) {
		Box bx = loc_it->second;
		bx.widen();
		result.insert(make_pair<DiscreteLocation,Box>(loc_it->first,bx));
	}

	return result;
}


HybridBoxes
unbounded_hybrid_boxes(const HybridSpace& hspace)
{
	HybridBoxes result;

	for (HybridSpace::const_iterator space_it = hspace.begin(); space_it != hspace.end(); ++space_it)
		result.insert(std::pair<DiscreteLocation,Box>(space_it->first,unbounded_box(space_it->second)));

	return result;
}


Box
project(const HybridBoxes& box, const std::vector<uint>& dimensions)
{
	ARIADNE_ASSERT_MSG(dimensions.size()>0, "Provide at least one dimension to project to.");

	Box result = Box::empty_box(dimensions.size());

	for (HybridBoxes::const_iterator loc_it = box.begin(); loc_it != box.end(); ++loc_it)
		result = hull(result,loc_it->second.project(dimensions));

	return result;
}


HybridBoxes
project(const Box& box, const std::vector<uint>& dimensions, const HybridSpace& target_space)
{
	ARIADNE_ASSERT_MSG(dimensions.size() == box.size(), "The original box and the projection sizes do not match.");

	HybridBoxes result = unbounded_hybrid_boxes(target_space);

	for (HybridBoxes::iterator loc_it = result.begin(); loc_it != result.end(); ++loc_it)
		for (uint i=0; i<dimensions.size();i++)
			loc_it->second[dimensions[i]] = box[i];

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

bool restricts(const HybridDenotableSet& hgts1, const HybridDenotableSet& hgts2) {
	ARIADNE_ASSERT_MSG(hgts1.grid() == hgts2.grid(), "For a restriction the grids of the two sets must be equal.");
    for(HybridDenotableSet::locations_const_iterator loc_iter=hgts1.locations_begin();
    		loc_iter!=hgts1.locations_end(); ++loc_iter) {
        if (restricts(loc_iter->second,hgts2.find(loc_iter->first)->second))
        	return true;
    }
    return false;
}


tribool disjoint(const HybridConstraintSet& cons_set, const HybridDenotableSet& grid_set)
{
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


HybridDenotableSet possibly_overlapping_subset(const HybridDenotableSet& hdenotable_set, const HybridConstraintSet& cons_set)
{
	HybridDenotableSet result(hdenotable_set.grid());

    for (HybridDenotableSet::locations_const_iterator gts_it = hdenotable_set.locations_begin(); gts_it != hdenotable_set.locations_end(); ++gts_it) {
    	HybridConstraintSet::locations_const_iterator cs_it = cons_set.find(gts_it->first);
    	if (cs_it != cons_set.locations_end())
    		result[gts_it->first] = possibly_overlapping_subset(gts_it->second,cs_it->second);
    	else
    		result[gts_it->first] = gts_it->second;
    }

	return result;
}


HybridDenotableSet definitely_covered_subset(const HybridDenotableSet& hdenotable_set, const HybridConstraintSet& cons_set)
{
	HybridDenotableSet result(hdenotable_set.grid());

    for (HybridDenotableSet::locations_const_iterator gts_it = hdenotable_set.locations_begin(); gts_it != hdenotable_set.locations_end(); ++gts_it) {
    	HybridConstraintSet::locations_const_iterator cs_it = cons_set.find(gts_it->first);
    	if (cs_it != cons_set.locations_end())
    		result[gts_it->first] = definitely_covered_subset(gts_it->second,cs_it->second);
    	else
    		result[gts_it->first] = gts_it->second;
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
			gts_it != grid_set.locations_end(); ++gts_it) {
		result.adjoin(project_down(gts_it->second,indices));
	}

	return result;
}


} // namespace Ariadne
