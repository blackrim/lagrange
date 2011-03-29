/*
 * node.h
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 */

#ifndef NODE_H_
#define NODE_H_

#include <string>
#include <vector>
#include <map>
using namespace std;


#include "node_object.h"
#include "vector_node_object.h"
#include "BranchSegment.h"

class Node{
private:
	double BL;//branch lengths
	double height; // could be from tip or from root
	int number;
	string name;
	Node * parent;
	vector<Node *> children;
	map<string,NodeObject *> assoc;
	map<string, vector<double> > assocDV;
	string comment;
	vector<BranchSegment> * segs;

public:
	Node();
	Node(Node * parent);
	Node(double bl,int number,string name, Node * parent);

	vector<Node*> getChildren();
	bool isExternal();
	bool isInternal();
	bool isRoot();
	bool hasParent();
	void setParent(Node & p);
	int getNumber();
	void setNumber(int n);
	double getBL();
	void setBL(double bl);
	double getHeight();
	void setHeight(double he);
	bool hasChild(Node & test);
	bool addChild(Node & c);
	bool removeChild(Node & c);
	Node & getChild(int c);
	string getName();
	string getComment();
	void setName(string s);
	void setComment(string s);
	string getNewick(bool bl);
	string getNewickOBL(string obj);
	string getNewick(bool bl,string obj);
	Node * getParent();
	int getChildCount();
	void assocObject(string name,NodeObject & obj);
	void assocDoubleVector(string name, vector<double> & obj);
	vector<double> * getDoubleVector(string name);
	void deleteDoubleVector(string name);
	void assocSegVector(vector<BranchSegment> & obj);
	vector<BranchSegment> * getSegVector();
	void deleteSegVector();
	NodeObject * getObject(string name);

};

#endif /* NODE_H_ */
