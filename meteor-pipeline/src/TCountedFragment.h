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

#ifndef TCOUNTEDFRAGMENT_H_
#define TCOUNTEDFRAGMENT_H_

#include "TSAMDefinition.h"

  const unsigned char C_NO_READ = 0; // 0000
  const unsigned char C_MULTIPLE_UNDIRECT_READ = 1; // 0001
  const unsigned char C_MULTIPLE_DIRECT_READ = 2; // 0010
  const unsigned char C_UNIQUE_UNDIRECT_READ = 4; // 0100
  const unsigned char C_UNIQUE_DIRECT_READ = 8; // 1000

  typedef std::vector<unsigned char> TByteArray; //// was array of byte (delphi)

  //// not used
//  TFragmentProfileStatsItem = record
//    long long TotalReadCount;
//    long long UniqueReadCount;
//    double SharedReadCount;
//    double SmartSharedReadCount;
//  end;

  //// not used
//  TFragmentProfileItem = record
//    /* Allow to see the coverage of the fragment according read type by bit-bit operations (OR)
//      C_NO_READ: byte = 0;  // 0000
//      C_MULTIPLE_UNDIRECT_READ: byte = 1; // 0001
//      C_MULTIPLE_DIRECT_READ: byte = 2; // 0010
//      C_UNIQUE_UNDIRECT_READ: byte = 4; // 0100
//      C_UNIQUE_DIRECT_READ: byte = 8; // 1000
//    */
//    byte Status;
//
//    TFragmentProfileStatsItem BothStrand;
//    TFragmentProfileStatsItem ForwardStrand;
//    TFragmentProfileStatsItem ReverseStrand;
//  end;

  class TCountedFragment {
  private:
    long long FFragmentID;
    long long FFragmentStartLocation;
    long long FFragmentEndLocation;
    long long FFragmentSize;
    bool FFragmentStrand;
    long long FTotalReadCount;
    double FSharedReadCount;
    double FSmartSharedReadCount;
    long long FUniqueReadCount;
    long long FDirectTotalReadCount;
    double FDirectSharedReadCount;
    double FDirectSmartSharedReadCount;
    long long FDirectUniqueReadCount;

    long long FTotalReadSumCoverage;
    double FSharedReadSumCoverage;
    double FSmartSharedReadSumCoverage;
    long long FUniqueReadSumCoverage;
    long long FDirectTotalReadSumCoverage;
    double FDirectSharedReadSumCoverage;
    double FDirectSmartSharedReadSumCoverage;
    long long FDirectUniqueReadSumCoverage;

    double FTotalReadMeanCoverage;
    double FSharedReadMeanCoverage;
    double FSmartSharedReadMeanCoverage;
    double FUniqueReadMeanCoverage;
    double FDirectTotalReadMeanCoverage;
    double FDirectSharedReadMeanCoverage;
    double FDirectSmartSharedReadMeanCoverage;
    double FDirectUniqueReadMeanCoverage;

    long long FUnDirectTotalReadSumCoverage;
    double FUnDirectSharedReadSumCoverage;
    double FUnDirectSmartSharedReadSumCoverage;
    long long FUnDirectUniqueReadSumCoverage;
    double FUnDirectTotalReadMeanCoverage;
    double FUnDirectSharedReadMeanCoverage;
    double FUnDirectSmartSharedReadMeanCoverage;
    double FUnDirectUniqueReadMeanCoverage;

    double FTotalReadCoverage;
    double FUniqueReadCoverage;
    double FDirectTotalReadCoverage;
    double FDirectUniqueReadCoverage;
    double FUnDirectTotalReadCoverage;
    double FUnDirectUniqueReadCoverage;

    TByteArray FFragmentProfile;
    /* Allow to see the coverage of the fragment according read type by bit-bit operations (OR)
      C_NO_READ: byte = 0;  // 0000
      C_MULTIPLE_UNDIRECT_READ: byte = 1; // 0001
      C_MULTIPLE_DIRECT_READ: byte = 2; // 0010
      C_UNIQUE_UNDIRECT_READ: byte = 4; // 0100
      C_UNIQUE_DIRECT_READ: byte = 8; // 1000
    */

  public:
    long long GetFragmentID() const;
    long long GetFragmentSize() const;
    long long GetFragmentStartLocation() const;
    long long GetFragmentEndLocation() const;
    long long GetTotalReadCount() const;
    double GetSharedReadCount() const;
    double GetSmartSharedReadCount() const;
    long long GetUniqueReadCount() const;
    long long GetDirectTotalReadCount() const;
    double GetDirectSharedReadCount() const;
    double GetDirectSmartSharedReadCount() const;
    long long GetDirectUniqueReadCount() const;

    long long GetTotalReadSumCoverage() const;
    double GetSharedReadSumCoverage() const;
    double GetSmartSharedReadSumCoverage() const;
    long long GetUniqueReadSumCoverage() const;
    long long GetDirectTotalReadSumCoverage() const;
    double GetDirectSharedReadSumCoverage() const;
    double GetDirectSmartSharedReadSumCoverage() const;
    long long GetDirectUniqueReadSumCoverage() const;

    double GetTotalReadCoverage() const;
    double GetUniqueReadCoverage() const;
    double GetDirectTotalReadCoverage() const;
    double GetDirectUniqueReadCoverage() const;
    double GetUnDirectTotalReadCoverage() const;
    double GetUnDirectUniqueReadCoverage() const;

    double GetTotalReadMeanCoverage() const;
    double GetSharedReadMeanCoverage() const;
    double GetSmartSharedReadMeanCoverage() const;
    double GetUniqueReadMeanCoverage() const;
    double GetDirectTotalReadMeanCoverage() const;
    double GetDirectSharedReadMeanCoverage() const;
    double GetDirectSmartSharedReadMeanCoverage() const;
    double GetDirectUniqueReadMeanCoverage() const;

    long long GetUnDirectTotalReadSumCoverage() const;
    double GetUnDirectSharedReadSumCoverage() const;
    double GetUnDirectSmartSharedReadSumCoverage() const;
    long long GetUnDirectUniqueReadSumCoverage() const;

    double GetUnDirectTotalReadMeanCoverage() const;
    double GetUnDirectSharedReadMeanCoverage() const;
    double GetUnDirectSmartSharedReadMeanCoverage() const;
    double GetUnDirectUniqueReadMeanCoverage() const;

    long long GetUnDirectTotalReadCount() const;
    double GetUnDirectSharedReadCount() const;
    double GetUnDirectSmartSharedReadCount() const;
    long long GetUnDirectUniqueReadCount() const;

    const TByteArray & GetFragmentProfile() const;

    TCountedFragment(
		const long long aFragmentID,
		const long long aFragmentStartLocation,
		const long long aFragmentEndLocation,
		const long long aFragmentSize,
      const unsigned char aFragmentStrand); //overload;

    TCountedFragment(); //overload;
    virtual ~TCountedFragment(); //override;

    void AddRead(
		const long long aReadStartLocation,
		const long long aReadEndLocation,
		const bool aReadStrand,
		const long long aReadMatchesCount); //overload;
    // aCigarArray describes the alignment when SAM output
    void AddRead(
		const int aReadStartLocation,
		const int aReadEndLocation,
		const bool aReadStrand,
		const int aReadMatchesCount,
		const TCigarTypeArray & aCigarArray); //overload;

    void Finalize();

    /*
      Two next procedures allow to calculate smart distribution of multiple reads
      taking account unique read fraction
    */
    void UpdateSmartSharedCount(const long long aMappedReadLength, const bool aReadStrand, const double aSmartSharedCount); //overload;
    // aCigarArray describes the alignment when SAM output
    void UpdateSmartSharedCount(const TCigarTypeArray & aCigarArray, const bool aReadStrand, const double aSmartSharedCount); //overload;
};


#endif /* TCOUNTEDFRAGMENT_H_ */
