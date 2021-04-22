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

#ifndef TPROFILEINIFILE_H_
#define TPROFILEINIFILE_H_

#include "MeteorConstant.h"
#include "IniFileManager.h"

class TProfileIniFile: public IniFileManager {

	private:
    	bool FReadOnly;

	public:
		TProfileIniFile(const std::string & filename, bool aIsReadOnly = false);
		virtual ~TProfileIniFile();

		//void SetProfileDate(const TDateTime & Value);
		void SetProfileDate(const std::string & Value);
		void SetReferenceName(const std::string & Value);
		void SetProjectName(const std::string & Value);
		void SetCountedReadType(const TCountedReadType & Value);
		void SetCountedReadType(const std::string & Value);
		void SetNormalizationType(const TMeteorProfileNormalizationType & Value);
		void SetDatabaseType(const TMeteorDBType & Value);

		void SetProfileFileName(const std::string & Value);
		void SetProfileName(const std::string & Value);
		void SetMeteorConfigPath(const std::string & Value);

		std::string GetConditionName() const;
		void SetConditionName(const std::string & Value);

		void SetReadCleaningMethod(const TReadCleaningMethodType & Value);

		void SetReadCleaningParameters(const std::string & Value);

		std::string GetSampleName() const;
		void SetSampleName(const std::string & Value);

		std::string GetProfiledSampleFullName(int Index) const;

		void SetProfileDirectory(const std::string & Value);

		int AddProfiledSample(const std::map<std::string, std::string> & aMapSection);
		int AddProfiledSample(const std::string & aSampleName, const std::string & aSampleFullName, const std::string & aConditionName, const std::string & aSequencingDate);
		//// TODO if needed: overloaded versions of AddProfiledSample() if needed

		void SetCountingStatistics(const int aProfiledSampleIndex, const TCountingStatistics & aCountingStatistics);
		TCountingStatistics GetCountingStatistics(const int aProfiledSampleIndex) const;

		TCountingStatistics GetCountingStatistics() const;
		void SetCountingStatistics(const TCountingStatistics & aCountingStatistics);

		void SetExcludedReferenceCountingStatistics(const TStrings & aExcludedReferenceNameList, const std::vector<int> & aExcludedReferenceReadCountArray);
		void SetExcludedReferenceCountingStatistics(const int aProfiledSampleIndex, const TStrings & aExcludedReferenceNameList, const std::vector<int> & aExcludedReferenceReadCountArray);

		void SetExcludedReferenceNameList(const TStrings & aExcludedReferenceNameList);

		void SetIndexedSequencedBaseCount(const int64 Value);
		void SetIndexedSequencedReadCount(const int Value);

		void SetIntersectCount(const int aProfiledSampleIndex, const int aIntersectDescrIndex, const int aCountValue);
		void SetIntersectCount(const int aIntersectDescrIndex, const int aCountValue);


		bool GetIsBestAlignment() const;
		void SetIsBestAlignment(const bool Value);

		void SetIsMappingMismatchesPercentage(const bool Value);
		void SetMappedReadLength(const int Value);
		void SetMappedReadLengthType(const TMappedReadLengthType Value);
		void SetMappingCmdLine(const std::string & Value);
		void SetMappingMatchesCount(const int Value);
		void SetMappingMismatchesCount(const int Value);

		int GetSoftClippingLengthCutoff() const; //// LOCAL
		int GetAlignmentLengthCutoff() const; //// LOCAL
	    bool GetIsLocalMapping() const;
	    bool GetKeepInternalLocalAlignment() const;
	    void SetSoftClippingLengthCutoff(const int Value);
	    void SetAlignmentLengthCutoff(const int Value);
	    void SetIsLocalMapping(const bool Value);
	    void SetKeepInternalLocalAlignment(const bool Value);

		std::string GetMappingToolName() const;
		void SetMappingToolName(const std::string & Value);

		std::string GetProjectName() const;
		TReadCleaningMethodType GetReadCleaningMethod() const;
		std::string GetReadCleaningParameters() const;

		int GetSequencedRawReadCount() const;
		void SetSequencedRawReadCount(const int Value);

		int GetSequencedHostReadCount() const;
		void SetSequencedHostReadCount(const int Value);

		int GetSequencedReadCount() const;
		void SetSequencedReadCount(const int Value);

		int GetProfiledSampleCount() const;

		bool GetProfiledSample(const int aProfileSampleIndex, std::map<std::string, std::string> & aMapSection) const;

		std::string GetReferenceName() const;
		int GetMappingMatchesCount() const;
		int GetMappingMismatchesCount() const;
		int GetMappedReadLength() const;
		bool GetIsMappingMismatchesPercentage() const;
		std::string GetProfiledSequencingDate(int Index) const;
		TMappedReadLengthType GetMappedReadLengthType() const;
		std::string GetMappingCmdLine() const;
		void GetExcludedReferenceNameList(TStrings & aExcludedReferenceNameList) const;

		int GetExcludedReferenceReadCount(const int aExcludedReferenceIndex) const;
		int GetExcludedReferenceReadCount(const int aProfiledSampleIndex, const int aExcludedReferenceIndex) const;
};

#endif /* TPROFILEINIFILE_H_ */
