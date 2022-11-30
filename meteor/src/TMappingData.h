/*
 * Copyright 2017-2020, Franck Gauthier <franck.gauthier@inrae.fr>, Nicolas Pons <nicolas.pons@inrae.fr>
 *
 * This file is part of Meteor v3.2.
 *
 * Meteor v3.2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Meteor v3.2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Meteor v3.2. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef TMAPPINGDATA_H_
#define TMAPPINGDATA_H_

#include <string>
#include <vector>
#include "MeteorConstant.h"
#include "TCountedFragmentList.h"
#include "TReferenceAnnotationDatabase.h"

//using namespace std;

typedef unsigned char TIntersectionStatus; // => 7 intersected reference max (7+1 bits)

const TIntersectionStatus C_MAIN_MAPPING_INTERSECTION_STATUS = 1;
const TIntersectionStatus C_NO_MAPPING_INTERSECTION_STATUS = 0;

typedef struct  {
	unsigned short MatchesCount;
	double MinDistance; // used to get the best alignment only.
	bool IsRejected;
	TIntersectionStatus IntersectionStatus;
} TMappingIntersectionItem;

typedef std::vector<TMappingIntersectionItem> TMappingIntersectionArray;

class TMappingData {

	private:
		void CreateMappingIntersectionArray(const int aMappedReadCount);

	protected:
	    int FMappedReadCount;
	    double FMaxDistance;
	    bool FIsRelativeDistance;
	    bool FIsLocalAlignment;
	    bool FKeepInternalLocalAlignment;
	    int FAlignmentLengthCutoff;
	    int FSoftClippingLengthCutoff;
	    TMappingFileFormat FMappingFileFormat;
	    TMappingIntersectionArray FMappingIntersectionArray;
	    std::vector<int> FExcludedReadCountArray;
	    int FExcludedReadCountArraySize;
	    TReferenceAnnotationDatabase * FReferenceAnnotationDatabase; //// not owned, do not delete in destructor.

	    void AddSAMExcludedMapping(
	    		const std::string & aSAMFileName, const int & aIndexExcludedMapping,
				const double & aExcludedMaxDistance, const bool aIsRelativeMaxDistance, const bool aIsOnBothStrand); ////
	    void AddBowtieExcludedMapping(
	    		const std::string & aFileName, const int & aIndexExcludedMapping,
				const double & aExcludedMaxDistance, const bool aIsRelativeMaxDistance, const bool aIsOnBothStrand); ////
	public:

//		TMappingData();

	    TMappingData(
	    		const int aMappedReadCount, const double aMaxDistance, const bool aIsRelativeDistance,
	    		const bool aIsLocalAlignment, const bool aKeepInternalLocalAlignment,
	    		const int aAlignmentLengthCutoff, const int aSoftClippingLengthCutoff, const TMappingFileFormat & aMappingFileFormat);
		virtual ~TMappingData();

		virtual void AddMainMapping(const std::string & aFileName) = 0; // pur virtual method
		void AddExcludedMapping(
				const std::string & aFileName,
				const int aIndexExcludedMapping, const double aExcludedMaxDistance,
				const bool aIsRelativeMaxDistance, const bool aIsOnBothStrand);

		virtual void CountFragment(TCountedFragmentList & aCountedFragmentList, TReferenceAnnotationDatabase & aReferenceAnnotationDatabase) = 0;
		virtual void UpdateSmartSharedCount(TCountedFragmentList & aCountedFragmentList, TReferenceAnnotationDatabase & aReferenceAnnotationDatabase) = 0; //// LOCAL
		void IncCountingStatistics(TCountingStatistics & aCountingStatistics) const;
		int GetExcludedReadCount(const int aIndexExcludedMapping) const;
		void IncIntersectionCount(std::vector<int> & aIntersectionCountArray);
};

#endif /* TMAPPINGDATA_H_ */
