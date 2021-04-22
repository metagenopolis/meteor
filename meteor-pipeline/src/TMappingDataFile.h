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

#ifndef TMAPPINGDATAFILE_H_
#define TMAPPINGDATAFILE_H_

#include <string>
#include "MeteorConstant.h"
#include "TCountedFragment.h"
#include "TCountedFragmentList.h"
#include "TReferenceAnnotationDatabase.h"
#include "TMappingData.h"
#include "TReadMatchInfoList.h"
#include "TSAMReader.h"
#include "TSAMDefinition.h"

//using namespace std;

// class managing mapping result directly with file

class TMappingDataFile: public TMappingData {

private:
	std::string FMainMappingFileName;

	void CountFragmentFromBowtieFile(const TCountedFragmentList & aCountedFragmentList);
	void CountFragmentFromSAMFile_old(TCountedFragmentList & aCountedFragmentList);
	void CountFragmentFromSAMFile(TCountedFragmentList & aCountedFragmentList);
	void SetMappingIntersectionArray(TReadMatchInfoList & aReadMatchInfoList);
	void CountSAMReadList(TReadMatchInfoList & aReadMatchInfoList, TCountedFragmentList & aCountedFragmentList);
	void UpdateSmartSharedCountFromBowtieFile(const TCountedFragmentList & aCountedFragmentList);
	void UpdateSmartSharedCountFromSAMFile(TCountedFragmentList & aCountedFragmentList, TReferenceAnnotationDatabase & aReferenceAnnotationDatabase); //// LOCAL
	void UpdateSmartSharedCountSAMReadList(TReadMatchInfoList & aReadMatchInfoList, TCountedFragmentList & aCountedFragmentList);
	void TagReadMatchInfoList(TReadMatchInfoList & aReadMatchInfoList);

public:
//	TMappingDataFile();

    TMappingDataFile(
    		const int aMappedReadCount, const double aMaxDistance, const bool aIsRelativeDistance,
    		const bool aIsLocalAlignment, const bool aKeepInternalLocalAlignment,
    		const int aAlignmentLengthCutoff, const int aSoftClippingLengthCutoff, const TMappingFileFormat & aMappingFileFormat);

	virtual ~TMappingDataFile();

	void AddMainMapping(const std::string & aFileName); // inherited //// LOCAL
	void CountFragment(TCountedFragmentList & aCountedFragmentList, TReferenceAnnotationDatabase & aReferenceAnnotationDatabase); // inherited
	void UpdateSmartSharedCount(TCountedFragmentList & aCountedFragmentList, TReferenceAnnotationDatabase & aReferenceAnnotationDatabase); //// LOCAL
};

#endif /* TMAPPINGDATAFILE_H_ */
