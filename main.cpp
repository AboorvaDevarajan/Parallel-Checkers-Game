#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <strings.h>
#include <vector>
#include <limits>
#include<mpi.h>

#include "Test.h"
#include "Player.h"

using namespace std;
int initial_depth = 3;
int inc_depth = 4;

/*
* Indicates the max and min values for tree search Alpha and beta (-/+ Infinity)
*/
int imin = std::numeric_limits<int>::min(); 
int imax = std::numeric_limits<int>::max(); 


/*
*Function to parse the file and read the test Tests
*/
void readTestTests(char *fileName, vector<Test*>*Tests) {
	string inputString;
	ifstream myfile(fileName);
	if (myfile.is_open()) {
		while (myfile.good() && (getline(myfile, inputString) != NULL)) {
			if (strncmp(inputString.c_str(), "Test", 4) == 0) {
				Test *newTest = new Test();

				for (int i = 0; i < 8; ++i) {
					if (myfile.good()
							&& (getline(myfile, inputString) != NULL)) {
						for (int j = 0; j < 8; ++j) {
							if ((inputString[j] == '+')
									|| (inputString[j] == 'O')) {
								newTest->checkerBoard[i][j] = 0;
							} else if (inputString[j] == 'A') {
								newTest->checkerBoard[i][j] = 1;
							} else if (inputString[j] == 'k') {
								newTest->checkerBoard[i][j] = 2;
							} else if (inputString[j] == 'B') {
								newTest->checkerBoard[i][j] = -1;
							} else if (inputString[j] == 'K') {
								newTest->checkerBoard[i][j] = -2;
							}
						}
					}
				}
				Tests->push_back(newTest);
			}
		}
	}
}

/*
*  Parsing the command line input. And check for correctness.
*/
char* parseCommadInput(int argc, char *argv[], vector<Test*> *Tests) {
	if (argc < 2) {
		cerr << "Invalid command. Correct format: proj_name input.txt output.txt"
				<< endl;
		exit(1);
	}
	
	readTestTests(argv[0], Tests);
	return argv[1];
}


/*
*   Generating Heuristic value for a given node.
*/
int getHeuristicValue(Test * newTest) {
	int h = 0;
	for (int i = 0; i < 8; ++i) {
		for (int j = 0; j < 8; ++j) {
			h += newTest->checkerBoard[i][j];
		}
	}
	return h;
}

/*
*   Function to copy array from one to another
*/
void copyArray(int arrayFrom[8][8], int arrayTo[8][8]) {
	for (int i = 0; i < 8; ++i) {
		for (int j = 0; j < 8; ++j) {
			arrayTo[i][j] = arrayFrom[i][j];
		}
	}
}

/*
*   function to print the children nodes
*/
void printChildren(vector<Test*> *children) {
	std::vector<Test*>::iterator itr;
	for (itr = children->begin(); itr < children->end(); ++itr) {
		for (int i = 0; i < 8; ++i) {
			for (int j = 0; j < 8; ++j) {
				cout << (*itr)->checkerBoard[i][j];
			}
			cout << endl;
		}
		cout << "------------------------" << endl;
	}
}

/*
*  print the tree search process
*/
void printAlphaBeta(int alpha, int beta, ofstream *output) {
	std::stringstream ss;
	ss << alpha;
	std::string alphaStr(ss.str());
	std::stringstream ss1;
	ss1 << beta;
	std::string betaStr(ss1.str());
	if (alpha == imin) {
		alphaStr = std::string("-Inf");
	} else if (alpha == imax) {
		alphaStr = std::string("Inf");
	}
	if (beta == imin) {
		betaStr = std::string("-Inf");
	} else if (beta == imax) {
		betaStr = std::string("Inf");
	}

	*output << "alpha=" << alphaStr << ", beta=" << betaStr << "." << endl;
}


/*
*  Finding all possible children for a given node.
*  Initially checks for jumps if any, which has the higher priority and then goes for the single move.
*/
void findAllPossibleMoves(Test* newTest, vector<Test*> *children, bool isMax) {
	//Find all the Tests of jumps if any
	//This is same as finding all the Tests where opponent (-1 or -2) is
	//on the diagonally opposite position. Spl Test for kings
	if (isMax) {
		for (int i = 0; i < 8; ++i) {
			for (int j = 0; j < 8; ++j) {
				if (newTest->checkerBoard[i][j] == 1) {
					//Pawns
					if ((i > 1) && (j < 6)
							&& ((newTest->checkerBoard[(i - 1)][(j + 1)] == -1)
									|| (newTest->checkerBoard[i - 1][j + 1] == -2))) // right jump check
						{
						if (newTest->checkerBoard[i - 2][j + 2] == 0) {
							Test *nBoard = new Test();
							copyArray(newTest->checkerBoard, nBoard->checkerBoard);
							nBoard->checkerBoard[i - 2][j + 2] =
									nBoard->checkerBoard[i][j];
							nBoard->checkerBoard[i - 1][j + 1] = 0;
							nBoard->checkerBoard[i][j] = 0;
							nBoard->mFrom=(i*10)+j;
							nBoard->mTo=((i - 2)*10)+(j + 2);	
							children->push_back(nBoard);
						}
					}
					if ((i > 1) && (j > 1)
							&& ((newTest->checkerBoard[i - 1][j - 1] == -1)
									|| (newTest->checkerBoard[i - 1][j - 1] == -2)))//left jump check
						 {
						if (newTest->checkerBoard[i - 2][j - 2] == 0) {
							Test *nBoard = new Test();
							copyArray(newTest->checkerBoard, nBoard->checkerBoard);
							nBoard->checkerBoard[i - 2][j - 2] =
									nBoard->checkerBoard[i][j];
							nBoard->checkerBoard[i - 1][j - 1] = 0;
							nBoard->checkerBoard[i][j] = 0;
							nBoard->mFrom=(i*10)+j;
							nBoard->mTo=((i-2)*10)+(j-2);
							children->push_back(nBoard);
						}
					}
				} else if (newTest->checkerBoard[i][j] == 2) {
					//Kings
					//rightside forward jump
					if ((i > 1) && (j < 6)
							&& ((newTest->checkerBoard[i - 1][j + 1] == -1)
									|| (newTest->checkerBoard[i - 1][j + 1]) == -2)) {
						if (newTest->checkerBoard[i - 2][j + 2] == 0) {
							Test *nBoard = new Test();
							copyArray(newTest->checkerBoard, nBoard->checkerBoard);
							nBoard->checkerBoard[i - 2][j + 2] =
									nBoard->checkerBoard[i][j];
							nBoard->checkerBoard[i - 1][j + 1] = 0;
							nBoard->checkerBoard[i][j] = 0;
							nBoard->mFrom=(i*10)+j;
							nBoard->mTo=((i-2)*10)+(j+2);
							children->push_back(nBoard);
						}
					}
					//leftside forward jump
					if ((i > 1) && (j > 1)
							&& ((newTest->checkerBoard[i - 1][j - 1] == -1)
									|| (newTest->checkerBoard[i - 1][j - 1] == -2))) {
						if (newTest->checkerBoard[i - 2][j - 2] == 0) {
							Test *nBoard = new Test();
							copyArray(newTest->checkerBoard, nBoard->checkerBoard);
							nBoard->checkerBoard[i - 2][j - 2] =
									nBoard->checkerBoard[i][j];
							nBoard->checkerBoard[i - 1][j - 1] = 0;
							nBoard->checkerBoard[i][j] = 0;
							nBoard->mFrom=(i*10)+j;
							nBoard->mTo=((i-2)*10)+(j-2);
							children->push_back(nBoard);
						}
					}
					//right side backward jump
					if ((i < 6) && (j < 6)
							&& ((newTest->checkerBoard[i + 1][j + 1] == -1)
									|| (newTest->checkerBoard[i + 1][j + 1] == -2))) {
						if (newTest->checkerBoard[i + 2][j + 2] == 0) {
							Test *nBoard = new Test();
							copyArray(newTest->checkerBoard, nBoard->checkerBoard);
							nBoard->checkerBoard[i + 2][j + 2] =
									nBoard->checkerBoard[i][j];
							nBoard->checkerBoard[i + 1][j + 1] = 0;
							nBoard->checkerBoard[i][j] = 0;
							nBoard->mFrom=(i*10)+j;
							nBoard->mTo=((i+2)*10)+(j+2);
							children->push_back(nBoard);
						}
					}
					//leftside backward jump
					if ((i < 6) && (j > 1)
							&& ((newTest->checkerBoard[i + 1][j - 1] == -1)
									|| (newTest->checkerBoard[i + 1][j - 1] == -2))) {
						if (newTest->checkerBoard[i + 2][j - 2] == 0) {
							Test *nBoard = new Test();
							copyArray(newTest->checkerBoard, nBoard->checkerBoard);
							nBoard->checkerBoard[i + 2][j - 2] =
									nBoard->checkerBoard[i][j];
							nBoard->checkerBoard[i + 1][j - 1] = 0;
							nBoard->checkerBoard[i][j] = 0;
							nBoard->mFrom=(i*10)+j;
							nBoard->mTo=((i+2)*10)+(j-2);
							children->push_back(nBoard);
						}
					}
				}
			}
		}
	} else {
		for (int i = 0; i < 8; ++i) {
			for (int j = 0; j < 8; ++j) {
				if (newTest->checkerBoard[i][j] == -1) {
					//Pawns opponent
					if ((i < 6) && (j < 6)
							&& ((newTest->checkerBoard[(i + 1)][(j + 1)] == 1)
									|| (newTest->checkerBoard[i + 1][j + 1] == 2))) {
						if (newTest->checkerBoard[i + 2][j + 2] == 0) {
							Test *nBoard = new Test();
							copyArray(newTest->checkerBoard, nBoard->checkerBoard);
							nBoard->checkerBoard[i + 2][j + 2] =
									nBoard->checkerBoard[i][j];
							nBoard->checkerBoard[i + 1][j + 1] = 0;
							nBoard->checkerBoard[i][j] = 0;
							nBoard->mFrom=(i*10)+j;
							nBoard->mTo=((i+2)*10)+(j+2);
							children->push_back(nBoard);
						}
					}
					if ((i < 6) && (j > 1)
							&& ((newTest->checkerBoard[i + 1][j - 1] == 1)
									|| (newTest->checkerBoard[i + 1][j - 1] == 2))) {
						if (newTest->checkerBoard[i + 2][j - 2] == 0) {
							Test *nBoard = new Test();
							copyArray(newTest->checkerBoard, nBoard->checkerBoard);
							nBoard->checkerBoard[i + 2][j - 2] =
									nBoard->checkerBoard[i][j];
							nBoard->checkerBoard[i + 1][j - 1] = 0;
							nBoard->checkerBoard[i][j] = 0;
							nBoard->mFrom=(i*10)+j;
							nBoard->mTo=((i+2)*10)+(j-2);
							children->push_back(nBoard);
						}
					}
				} else if (newTest->checkerBoard[i][j] == -2) {
					//Kings
					if ((i > 1) && (j < 6)
							&& ((newTest->checkerBoard[i - 1][j + 1] == 1)
									|| (newTest->checkerBoard[i - 1][j + 1]) == 2)) {
						if (newTest->checkerBoard[i - 2][j + 2] == 0) {
							Test *nBoard = new Test();
							copyArray(newTest->checkerBoard, nBoard->checkerBoard);
							nBoard->checkerBoard[i - 2][j + 2] =
									nBoard->checkerBoard[i][j];
							nBoard->checkerBoard[i - 1][j + 1] = 0;
							nBoard->checkerBoard[i][j] = 0;
							nBoard->mFrom=(i*10)+j;
							nBoard->mTo=((i-2)*10)+(j+2);
							children->push_back(nBoard);
						}
					}
					if ((i > 1) && (j > 1)
							&& ((newTest->checkerBoard[i - 1][j - 1] == 1)
									|| (newTest->checkerBoard[i - 1][j - 1] == 2))) {
						if (newTest->checkerBoard[i - 2][j - 2] == 0) {
							Test *nBoard = new Test();
							copyArray(newTest->checkerBoard, nBoard->checkerBoard);
							nBoard->checkerBoard[i - 2][j - 2] =
									nBoard->checkerBoard[i][j];
							nBoard->checkerBoard[i - 1][j - 1] = 0;
							nBoard->checkerBoard[i][j] = 0;
							nBoard->mFrom=(i*10)+j;
							nBoard->mTo=((i-2)*10)+(j-2);
							children->push_back(nBoard);
						}
					}
					if ((i < 6) && (j < 6)
							&& ((newTest->checkerBoard[i + 1][j + 1] == 1)
									|| (newTest->checkerBoard[i + 1][j + 1] == 2))) {
						if (newTest->checkerBoard[i + 2][j + 2] == 0) {
							Test *nBoard = new Test();
							copyArray(newTest->checkerBoard, nBoard->checkerBoard);
							nBoard->checkerBoard[i + 2][j + 2] =
									nBoard->checkerBoard[i][j];
							nBoard->checkerBoard[i + 1][j + 1] = 0;
							nBoard->checkerBoard[i][j] = 0;
							nBoard->mFrom=(i*10)+j;
							nBoard->mTo=((i+2)*10)+(j+2);
							children->push_back(nBoard);
						}
					}
					if ((i < 6) && (j > 1)
							&& ((newTest->checkerBoard[i + 1][j - 1] == 1)
									|| (newTest->checkerBoard[i + 1][j - 1] == 2))) {
						if (newTest->checkerBoard[i + 2][j - 2] == 0) {
							Test *nBoard = new Test();
							copyArray(newTest->checkerBoard, nBoard->checkerBoard);
							nBoard->checkerBoard[i + 2][j - 2] =
									nBoard->checkerBoard[i][j];
							nBoard->checkerBoard[i + 1][j - 1] = 0;
							nBoard->checkerBoard[i][j] = 0;
							nBoard->mFrom=(i*10)+j;
							nBoard->mTo=((i+2)*10)+(j-2);
							children->push_back(nBoard);
						}
					}
				}
			}
		}
	}

	//If there are no Tests where you can jump then find all
	//the other possible moves
	//Find all moves for normal pawns
	//Find all possible moves for Kings
	if (children->empty()) {
		for (int i = 0; i < 8; ++i) {
			for (int j = 0; j < 8; ++j) {
				if (isMax && newTest->checkerBoard[i][j] == 1) {
					//Pawns
					//right side
					if ((i > 0) && (j < 7)
							&& (newTest->checkerBoard[(i - 1)][(j + 1)] == 0)) {
						Test *nBoard = new Test();
						copyArray(newTest->checkerBoard, nBoard->checkerBoard);
						nBoard->checkerBoard[i - 1][j + 1] =
								nBoard->checkerBoard[i][j];
						nBoard->checkerBoard[i][j] = 0;
						nBoard->mFrom=(i*10)+j;
						nBoard->mTo=((i-1)*10)+(j+1);		
						children->push_back(nBoard);
					}
					//left side
					if ((i > 0) && (j > 0)
							&& (newTest->checkerBoard[i - 1][j - 1] == 0)) {
						Test *nBoard = new Test();
						copyArray(newTest->checkerBoard, nBoard->checkerBoard);
						nBoard->checkerBoard[i - 1][j - 1] =
								nBoard->checkerBoard[i][j];
						nBoard->checkerBoard[i][j] = 0;
						nBoard->mFrom=(i*10)+j;
						nBoard->mTo=((i-1)*10)+(j-1);
						children->push_back(nBoard);
					}
				} else if ((!isMax) && newTest->checkerBoard[i][j] == -1) {
					//Pawns opponent
					//right side
					if ((i < 7) && (j < 7)
							&& (newTest->checkerBoard[(i + 1)][(j + 1)] == 0)) {
						Test *nBoard = new Test();
						copyArray(newTest->checkerBoard, nBoard->checkerBoard);
						nBoard->checkerBoard[i + 1][j + 1] =
								nBoard->checkerBoard[i][j];
						nBoard->checkerBoard[i][j] = 0;
						nBoard->mFrom=(i*10)+j;
						nBoard->mTo=((i+1)*10)+(j+1);	
						children->push_back(nBoard);
					}
					//left side
					if ((i < 7) && (j > 0)
							&& (newTest->checkerBoard[i + 1][j - 1] == 0)) {
						Test *nBoard = new Test();
						copyArray(newTest->checkerBoard, nBoard->checkerBoard);
						nBoard->checkerBoard[i + 1][j - 1] =
								nBoard->checkerBoard[i][j];
						nBoard->checkerBoard[i][j] = 0;
						nBoard->mFrom=(i*10)+j;
						nBoard->mTo=((i+1)*10)+(j-1);	
						children->push_back(nBoard);
					}
				} else if ((isMax && newTest->checkerBoard[i][j] == 2)
						|| ((!isMax) && newTest->checkerBoard[i][j] == -2)) {
					//Kings
					if ((i > 0) && (j < 7)
							&& (newTest->checkerBoard[i - 1][j + 1] == 0)) {
						Test *nBoard = new Test();
						copyArray(newTest->checkerBoard, nBoard->checkerBoard);
						nBoard->checkerBoard[i - 1][j + 1] =
								nBoard->checkerBoard[i][j];
						nBoard->checkerBoard[i][j] = 0;
						nBoard->mFrom=(i*10)+j;
						nBoard->mTo=((i-1)*10)+(j+1);	
						children->push_back(nBoard);
					}
					if ((i > 0) && (j > 0)
							&& (newTest->checkerBoard[i - 1][j - 1] == 0)) {
						Test *nBoard = new Test();
						copyArray(newTest->checkerBoard, nBoard->checkerBoard);
						nBoard->checkerBoard[i - 1][j - 1] =
								nBoard->checkerBoard[i][j];
						nBoard->checkerBoard[i][j] = 0;
						nBoard->mFrom=(i*10)+j;
						nBoard->mTo=((i-1)*10)+(j-1);	
						children->push_back(nBoard);
					}
					if ((i < 7) && (j < 7)
							&& (newTest->checkerBoard[i + 1][j + 1] == 0)) {
						Test *nBoard = new Test();
						copyArray(newTest->checkerBoard, nBoard->checkerBoard);
						nBoard->checkerBoard[i + 1][j + 1] =
								nBoard->checkerBoard[i][j];
						nBoard->checkerBoard[i][j] = 0;
						nBoard->mFrom=(i*10)+j;
						nBoard->mTo=((i+1)*10)+(j+1);	
						children->push_back(nBoard);
					}
					if ((i < 7) && (j > 0)
							&& (newTest->checkerBoard[i + 1][j - 1] == 0)) {
						Test *nBoard = new Test();
						copyArray(newTest->checkerBoard, nBoard->checkerBoard);
						nBoard->checkerBoard[i + 1][j - 1] =
								nBoard->checkerBoard[i][j];
						nBoard->checkerBoard[i][j] = 0;
						nBoard->mFrom=(i*10)+j;//nBoard->move << "(" << i << ", " << j << ") => ("
						nBoard->mTo=((i+1)*10)+(j-1);	//		<< (i + 1) << ", " << (j - 1) << ")";
						children->push_back(nBoard);
					}
				}
			}
		}
	}
}

/*
* alpha beta pruning algorithm
*/
 
int alphaBeta(Test *newTest, int depth, int alpha, int beta, Player *newPlayer,
		Player *notPlayer, ofstream *output) {
	if (depth == 0) {
		int temp = getHeuristicValue(newTest);
		newTest->value = temp;
		*output << "h=" << temp << "." << endl;
		return temp;
	}
	if (newPlayer->isMax) {
		vector<Test*> children;
		findAllPossibleMoves(newTest, &children, newPlayer->isMax);
		std::vector<Test*>::iterator itr;
		for (itr = children.begin(); itr < children.end(); ++itr) {
			for (int p = 0; p < initial_depth - depth; ++p) {
				*output << "|-";
			}
			*output << "A" << inc_depth - depth << ": " << (*itr)->mFrom<<"  to  "<<(*itr)->mTo << "."
					<< endl;
			alpha = max(alpha,
					alphaBeta(*itr, depth - 1, alpha, beta, notPlayer,
							newPlayer, output));
			if (beta <= alpha) {         //Beta cut-off
				++itr;
				if (itr < children.end()) {
					for (int p = 0; p < initial_depth - depth; ++p) {
						*output << "|-";
					}
					*output << "A" << inc_depth - depth << ": pruning ";
					//for each remaining children the branch is pruned
					int counter = 1;
					for (; itr < children.end(); ++itr) {
						if (counter == 1) {
							*output << (*itr)->mFrom<<"  to  "<<(*itr)->mTo ;
						} else {
							*output << ", " << (*itr)->mFrom<<"  to  "<<(*itr)->mTo ;
						}
						counter++;
					}
					*output << "; ";
					printAlphaBeta(alpha, beta, output);
				}
				break;
			}
		}
		newTest->value = imin;
		if (children.size() == 0) {
			newTest->value = getHeuristicValue(newTest);
			alpha = newTest->value;
		}
		for (itr = children.begin(); itr < children.end(); ++itr) {
			if ((*itr)->value > newTest->value) {
				newTest->value = (*itr)->value;
				newTest->bFrom=(*itr)->mFrom;
				newTest->bTo=(*itr)->mTo;
			}
		}
		if (depth == initial_depth) {
		cout<<"\n\n first move"<<newTest->mFrom<<"  to  "<<newTest->mTo 
					<< "." << endl;
			*output << endl << "first move: " << newTest->mFrom<<"  to  "<<newTest->mTo 
					<< "." << endl;
		}
		return alpha;
	} else {
		vector<Test*> children;
		findAllPossibleMoves(newTest, &children, newPlayer->isMax);
		std::vector<Test*>::iterator itr;
		for (itr = children.begin(); itr < children.end(); ++itr) {
			for (int p = 0; p < initial_depth - depth; ++p) {
				*output << "|-";
			}
			if (depth == 1) {
				*output << "B" << inc_depth - depth << ": " << (*itr)->mFrom<<"  to  "<<(*itr)->mTo 
						<< ";  ";
			} else {
				*output << "B" << inc_depth - depth << ": " << (*itr)->mFrom<<"  to  "<<(*itr)->mTo << "."
						<< endl;
			}
			beta = min(beta,
					alphaBeta(*itr, depth - 1, alpha, beta, notPlayer,
							newPlayer, output));
			if (beta <= alpha) {         //alpha cut-off
				++itr;
				if (itr < children.end()) {
					for (int p = 0; p < initial_depth - depth; ++p) {
						*output << "|-";
					}
					*output << "B" << inc_depth - depth << ": pruning ";
					//for each remaining children the branch is pruned
					int counter = 1;
					for (; itr < children.end(); ++itr) {
						if (counter == 1) {
							*output << (*itr)->mFrom<<"  to  "<<(*itr)->mTo ;
						} else {
							*output << ", " << (*itr)->mFrom<<"  to  "<<(*itr)->mTo ;
						}
						counter++;
					}
					*output << "; ";
					printAlphaBeta(alpha, beta, output);
				}
				break;
			}
		}
		newTest->value = imax;
		if (children.size() == 0) {
			newTest->value = getHeuristicValue(newTest);
			beta = newTest->value;
		}
		for (itr = children.begin(); itr < children.end(); ++itr) {
			if ((*itr)->value < newTest->value) {
				newTest->value = (*itr)->value;
				newTest->bFrom=(*itr)->mFrom;
				newTest->bTo=(*itr)->mTo;
			}
		}
		return beta;
	}
	return 0;
}

int main(int argc, char **argv) {
     int rank,size;
	MPI::Init(argc,argv);
	rank=MPI::COMM_WORLD.Get_rank();
	size=MPI::COMM_WORLD.Get_size();
	int count=0;
	int children_size;
	/*
	*   Parallelizing the Game tree using MPI functions
	*/
	
	if(rank==0){
	     vector<Test*> Tests;
	     char *out_name[2];
	     out_name[0]="input.txt";
	     out_name[1]="output.txt";
	
	     char* filename = parseCommadInput(3, out_name, &Tests);
	
	
	     std::vector<Test*>::iterator itr;
	     vector<Test*> children;
	     itr=Tests.begin();
	     findAllPossibleMoves(*itr,&children,1);
	
	     children_size=children.size();
		
	     int child_count;
	     int proc_to_assign=1;
          int tag=1;
          for(int send_child_count=1;send_child_count<size;send_child_count++)
          MPI::COMM_WORLD.Send(&children_size,sizeof(int),MPI::INT,send_child_count,111);
      
          std::vector<Test*>::iterator itr1;

		itr1=children.begin();

          cout<<"\n\n Children size  : "<<children_size;
          fflush(stdout);
      
      
          Test test;
      
          for(itr1=children.begin();itr1<children.end();itr1++){      
               test=**itr1;

     	     cout<<"\n\nAssigned tag values  "<<tag<<" and its processor: "<<proc_to_assign;
	          fflush(stdout);
	          MPI::COMM_WORLD.Send(&test,sizeof(Test),MPI::BYTE,proc_to_assign,tag);
	          if(proc_to_assign!=(size-1))
	               proc_to_assign++;
	          else{
	               proc_to_assign=1;
	          }
	          tag++;
          }
	}
	MPI::COMM_WORLD.Barrier();
	if(rank != 0){
          int assigned_child;
	     int rem;
	     int child_s;
          MPI::COMM_WORLD.Recv(&child_s,sizeof(int),MPI::INT,0,111);    
          assigned_child=child_s/(size-1);
          rem=child_s%(size-1);
          fflush(stdout);
          if(rem >= (rank)){
               assigned_child=assigned_child+1;
          }
          fflush(stdout);
          int tag_val;
          int child_count;
          Test sample[100];   //here occurs a conflict when two or more processes are allocated in the same process      
          tag_val=rank;
          cout<<"\nRank : "<<rank;
          fflush(stdout);
	     int tag=1000+rank;
          for(child_count=1;child_count<=assigned_child;child_count++){
         
               MPI::COMM_WORLD.Recv(&sample[child_count],sizeof(Test),MPI::BYTE,0,tag_val);
		     vector<Test*> Tests;
               Tests.push_back(&sample[child_count]);	
               ofstream output;
               output.open("output.txt");
               std::vector<Test*>::iterator itr;
               Test final_Test;
	
	          for (itr = Tests.begin(); itr < Tests.end(); ++itr) {
		          Player *A = new Player("A", true);
		          Player *B = new Player("B", false);
		          alphaBeta(*itr,initial_depth , imin, imax, B, A, &output);
		          final_Test=(**itr);
		          MPI::COMM_WORLD.Send(&final_Test,sizeof(Test),MPI::BYTE,0,tag);		
		          cout<<"\nafter alpha beta move value\n"<<"  "<<(*itr)->mFrom<<"   to  "<<(*itr)->mTo<<"   "<<(*itr)->value;
		
		          fflush(stdout);
		          //tag=tag+size-1;
	          }
	          tag=tag+size-1;
	          output.close();
                   
               cout<< "\n\nFrom (X,Y) : " << sample[child_count].mFrom ;
               fflush(stdout);
               cout<< "\n\nTo (X,Y) : " << sample[child_count].mTo ;
               fflush(stdout);
               tag_val=tag_val+(size-1);
               cout<<"\ntag value: "<<tag_val;
               fflush(stdout);
         
          } 
     } 
	MPI::COMM_WORLD.Barrier();
	if(rank==0)
	{
		int child_count;
		int proc_to_assign=1;
		int tag=1001;
		Test test[50];
		int itr_count;
		int test_count=0;
		for(itr_count=1;itr_count<size;itr_count++)
		{    tag=1000+itr_count;  
			int assigned_child=1;
			int rem;
			
			assigned_child=children_size/(size-1);
			rem=children_size%(size-1);
			
			if(rem >= (itr_count))
			{
				
				assigned_child=assigned_child+1;
			}
			
			 for(child_count=1;child_count<=assigned_child;child_count++)
			{
				
				MPI::COMM_WORLD.Recv(&test[test_count],sizeof(Test),MPI::BYTE,itr_count,tag);
				
				
				cout<<"\nIn rank 0"<< test[test_count].mFrom<<"  to  "<<test[test_count].mTo<<"  val  "<<test[test_count].value;
				
				tag=tag+(size-1);
				test_count++;
				
			}
			
			
          }
	     int h_count=0;
	     Test temp=test[h_count];
	     for(h_count=0;h_count<test_count;h_count++)
	     {
		     if(test[h_count].value>temp.value)
		     {
			     temp=test[h_count];
		     }
	     }
	     cout<<"\n\n Best move  "<<temp.mFrom<<"  to  "<<temp.mTo;
	
	}	
	
	MPI::Finalize();
	return 0;
}


