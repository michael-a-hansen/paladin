/*
 * The MIT License
 *
 * Copyright (c) 2016, 2017 Michael A. Hansen
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#ifndef SRC_MPI_UTIL_HPP_
#define SRC_MPI_UTIL_HPP_

#include <mpi.h>

namespace paladin {

struct MpiComm {
  int numRanks, myRank, rootRank;
  bool amRoot;

  MpiComm(int& argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    rootRank = 0;
    amRoot = (myRank == rootRank);
  }

  ~MpiComm() { MPI_Finalize(); }

  void barrier() const { MPI_Barrier(MPI_COMM_WORLD); }
};

/**
 * @brief convenience for MPI broadcasting a std::string
 * @param str a std::string to be broadcasted to slave ranks
 * @param myRank the current rank
 * @param numRanks the number of ranks
 * @param rootRank (optional) the index of the rootRank = 0 by default
 */
void broadcast_string(std::string& str, const MpiComm& comm) {
  if (comm.myRank == comm.rootRank) {
    for (int sendToRank = 0; sendToRank < comm.numRanks; ++sendToRank) {
      if (sendToRank != comm.rootRank) {
        MPI::COMM_WORLD.Send(str.c_str(), str.length(), MPI::CHAR, sendToRank,
                             0);
      }
    }
  } else {
    MPI_Status status;
    int size;
    MPI_Probe(comm.rootRank, 0, MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, MPI_CHAR, &size);

    std::vector<char> tmp(size);
    MPI_Recv(tmp.data(), size, MPI_CHAR, comm.rootRank, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    str.assign(tmp.begin(), tmp.end());
  }
  //  if (comm.myRank == comm.rootRank) {
  //    for (int sendToRank = 0; sendToRank < comm.numRanks; ++sendToRank)
  //      if (sendToRank != comm.rootRank)
  //        MPI_Send(str.c_str(), str.size(), MPI_CHAR, sendToRank, 0,
  //        MPI_COMM_WORLD);
  //  } else {
  //    MPI_Status status;
  //    int size;
  //    MPI_Probe(comm.rootRank, 0, MPI_COMM_WORLD, &status);
  //    MPI_Get_count(&status, MPI_CHAR, &size);
  //
  //    std::vector<char> tmp(size);
  //    MPI_Recv(tmp.data(), size, MPI_CHAR, comm.rootRank, 0, MPI_COMM_WORLD,
  //             MPI_STATUS_IGNORE);
  //    str.assign(tmp.begin(), tmp.end());
  //  }
}

}  // namespace paladin

#endif /* SRC_MPI_UTIL_HPP_ */
