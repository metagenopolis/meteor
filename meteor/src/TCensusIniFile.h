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

#ifndef TCENSUSINIFILE_H_
#define TCENSUSINIFILE_H_

#include "MeteorConstant.h"
#include "IniFileManager.h"

//using namespace std;


typedef struct {
    bool Exists;
    std::string Tag;
    std::string Name;
    bool IsCompressed;
    TSequenceFileFormat SequenceFileFormat;
    std::string CSFastaFileName, QualFastaFileName, FastQFileName;
//    switch(SequenceFileFormat){ //: TSequenceFileFormat ////
//    	case sffCSFastaAndQual:
//        //(CSFastaFileName, QualFastaFileName: std::string[255]);
//    	 std::string CSFastaFileName, QualFastaFileName; break;
//    	case sffFastQ:
//        std::string FastQFileName; break;
//    }
} TSampleFileInfo;

class TCensusIniFile : public IniFileManager {

private:
    bool FReadOnly;

public:
	TCensusIniFile(const std::string & filename, bool aIsReadOnly = false);
	virtual ~TCensusIniFile();

	void GetMappingBowtieFileNameList(TStrings & aMappingBowtieFileNameList);

	int GetProcessedReadCount() const;
	void SetProcessedReadCount(const int & Value);

	bool GetReadOnly() const;
	void SetReadOnly(bool value);

	TMappingFileFormat GetMappingFileFormat() const;
	void SetMappingFileFormat(const TMappingFileFormat & Value);

	std::string GetSequencingDate() const;

	std::string GetConditionName() const;
	std::string GetProjectName() const;
	std::string GetSampleName() const;

	std::string GetSequencingDevice() const;
	std::string GetFullSampleName() const;
	int GetIndexedSequencedReadCount() const;
	int GetSequencedReadCount() const;
	int GetSequencedHostReadCount() const;
	int64 GetIndexedSequenceBaseCount() const;
	TReadCleaningMethodType GetReadCleaningMethod() const;
	std::string GetReadCleaningParameters() const;

	bool GetIsCompressedSampleFile() const;
	void SetIsCompressedSampleFile(const bool & Value);

	int GetSequencedRawReadCount() const;

	bool GetIsLocalMapping() const ;
	void SetIsLocalMapping(const bool Value);
	int GetAlignmentLengthCutoff() const;
	int GetSoftClippingLengthCutoff() const;
	bool GetKeepInternalLocalAlignment() const;
	void SetAlignmentLengthCutoff(const int Value);
	void SetSoftClippingLengthCutoff(const int Value);
	void SetKeepInternalLocalAlignment(const bool Value);

	TCountingStatistics GetCountingStatistics() const;
	void GetExcludedReferenceNameList(TStrings & aExcludedReferenceNameList) const;
	std::string GetMappingTool() const;
	TSequenceFileFormat GetSequenceFileFormat() const;
	TSampleFileInfo GetSampleFileInfo() const;

};

#endif /* TCENSUSINIFILE_H_ */
