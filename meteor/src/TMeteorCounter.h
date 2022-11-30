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

#ifndef TMETEORCOUNTER_H_
#define TMETEORCOUNTER_H_

/* -----------------------------------------------------------------------------
 Unit Name: MeteorCounter.h
 Author:    Nicolas Pons (nicolas.pons@jouy.inra.fr
 (C) Copyright 2014 INRA MetaGenoPolis US1367 (http://www.mgps.eu)

 Meteor is under a proprietary software license
 Source code for the version 1 has been deposited to the APP
 (French Agency for the protection of softwares)
 APP number : IDDN.FR.001.420008.000.R.P.2013.000.30000

 Meteor source code includes AdvantageDB 11.x dependancy (http://devzone.advantagedatabase.com/)
 Please refer to the AdvantageDB's licence for any case of use

 Meteor calls the external software 'BOWTIE' for the read mapping task
 Please refer to the Bowtie's license (http://bowtie-bio.sourceforge.net/index.shtml)

 For any case of use, modification and/or distribution, you have to contact authors and MGP

 Purpose:
 Class allowing the gene counting for an indexed mapping

 History:
 ----------------------------------------------------------------------------- */
#include <vector>
#include <string>

#include "TCensusIniFile.h"
//#include "MappingDataset.h"
//#include "ReferenceIniFile.h"
#include "TNGSCensusData.h"
#include "MeteorConstant.h"
#include "TCountedFragmentList.h"
#include "TMeteorJobIniFile.h"
#include "TReferenceAnnotationDatabase.h"
#include "TMappingDataFile.h"
#include "TMappingData.h"


//using namespace std;
/*
 u_MappingDataset, u_CensusIniFile,
 u_ReferenceIniFile, u_NGSCensusData, u_MeteorConstant, u_CountedFragmentList,
 u_MeteorJobIniFile, u_ReferenceAnnotationDatabase, u_MappingDataFile,
 u_MappingData;
 */

class TMeteorCounter {

	private:
		// TMappingData is the base class of child classes:
		//   TMappingDataset (not implemented), was used on MS Windows platform (Delphi code) to deal with several mapping files (because of a bowtie2-build issue with large file)
		//   TMappingDataFile, used if only one mapping file exist
		// => It is necessary to use an array of POINTERS of TMappingData to use polymorphism of object.
		// => each pointer has to be allocated with operator new
		// => each pointer has to be de-allocated with operator delete (e.g. in the destructor)
		std::vector<TMappingData*> FMappingDataList; //// std::vector<TMappingData*> was TObjectList. tableau de pointeurs dans le code delphi.
		std::vector<TCountedReadType> FCountedReadTypeList;
		////TCountedFragmentList FCountedFragmentList; //// TODO uncomment if problem
		TReferenceAnnotationDatabase FReferenceAnnotationDatabase;

		std::string FCountingPrefixName;
		double FMaxDistance;
		bool FIsRelativeDistance;

		void GetCountingStatistics(TCountingStatistics & aCountingStatistics) const;
		int LoadReferenceAnnotationDatabase(const std::string & aReferenceIniFileName);
		void GetIntersectionCountArray(std::vector<int> & aIntersectionCountArray, const int aReferenceCount);

	public:
		TMeteorCounter(const std::vector<TCountedReadType> & aCountedReadTypeList);
		virtual ~TMeteorCounter();

		//DigestAlgorithm: TDigestAlgorithm; //// TODO ????
		void AddMappingLibrary(const std::string & aTmpDir, const std::string & aMainMappingCensusIniFileName, const double aMaxDistance,
				const bool aIsRelativeDistance,
				const bool aIsLocalAlignment, const bool aKeepInternalLocalAlignment,
				const int aAlignmentLengthCutoff, const int aSoftClippingLengthCutoff,
				const TStrings & aExcludedMappingCensusIniFileNameList,
				const TMappingReferencePropertiesArray & aExcludedMappingReferencePropertiesArray);

		void ProcessCounting(const std::string & aCountingPrefixName, const std::string & aNGSMappingDirectory,
				const std::string & aReferenceIniFileName, const std::string & aProfileIniFileName,
				const TStrings & aLibraryCensusIniFileNameList, const TStrings & aExcludedReferenceNameList,
				const TMappingReferenceProperties & aMainMappingReferenceProperties,
				const bool aOKToCacheData, const TMeteorDBType & aMeteorDBType);
};

#endif /* TMETEORCOUNTER_H_ */

