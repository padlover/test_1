#include <iostream>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#define infinte 2147483646
using namespace std;
class Node
{
public:
	typedef pair<int, int> point;
	Node(){};
	Node(int x, int y, int flow)
	{
	   coordinate.first = x;
	   coordinate.second = y;
	   maxflow = flow;
	   currentflow = 0;
     NodeIndex = -1;
     updated = false;
     vistited = false;
	}
	point coordinate;
   	int maxflow;
   	int currentflow;
    int NodeIndex;
    int minDist;
    bool updated;
    bool vistited;
};

bool myStr2Int(const string&, int&);
int abs(int);
int min(int&, int&);

int main(int argc, char** argv)
{
   //ifstream ifs(argv[1]);
   ofstream ofs(argv[2]);
   ifstream ifs(argv[1]);
   if(!ifs.is_open()) exit(-1);

   typedef pair<int, int> point;
   //file parsing
   vector<Node> sourceList;
   vector<Node> terminalList;
   
   int lineNum;
   int flow;
   int x, y;
   string str;
   ifs >> lineNum;
   for(int i = 0; i < lineNum; i++)
   {
      ifs >> x >> y;
      ifs >> str;
      if(!isdigit(str[0]) && str[0] == '+') str = str.substr(1);
      myStr2Int(str, flow);
      Node node(x, y, flow);
      if(node.maxflow > 0) sourceList.push_back(node);
      else if(node.maxflow < 0) terminalList.push_back(node);
   }
   //run greedy
     //(1) construct distance table
     int distanceTable[terminalList.size()][sourceList.size()];
     int flowTable[terminalList.size()][sourceList.size()];
     for(int i = 0; i < terminalList.size(); i++)
     {
     	for(int j = 0; j < sourceList.size(); j++)
     		flowTable[i][j] = 0;
     }

     //making distance form
     for(int i = 0; i < terminalList.size(); i++)
     {
       	for(int j = 0; j < sourceList.size(); j++)
       	{
            distanceTable[i][j] = abs(sourceList[j].coordinate.first - terminalList[i].coordinate.first)+abs(sourceList[j].coordinate.second - terminalList[i].coordinate.second);
       	}
     }
     typedef pair<int, int> edgeIndex;
     vector<edgeIndex> edges;
     for(int j = 0; j < sourceList.size(); j++)
     {
     	edges.clear();
     	for(int i = 0; i < terminalList.size(); i++)
     	{
            edges.push_back(edgeIndex(distanceTable[i][j], i));
     	}
     	std::sort(edges.begin(), edges.end());
     	while(sourceList[j].currentflow != sourceList[j].maxflow)
     	{
     		for(int k = 0; k < edges.size(); k++)
     		{
     			int index = edges[k].second;
     		    int remainflow = abs(terminalList[index].maxflow + terminalList[index].currentflow);
     			if(remainflow != 0)
     			{
     				flowTable[index][j] = (sourceList[j].maxflow - sourceList[j].currentflow > remainflow)? remainflow: sourceList[j].maxflow - sourceList[j].currentflow;
     			    sourceList[j].currentflow = sourceList[j].currentflow + flowTable[index][j];
     			    terminalList[index].currentflow = terminalList[index].currentflow + flowTable[index][j];
     			}
     		}
     	}
     }

//while negative cycle detected

   //construct residue graph
   int negativeRes[terminalList.size()][sourceList.size()];
   int positiveRes[terminalList.size()][sourceList.size()];
   int posResDist[terminalList.size()][sourceList.size()];
   int negaResDist[terminalList.size()][sourceList.size()];



//find negative cycle
   bool doupdate = true;
   vector<int> trace;
   int updateNum = 0;
   while(doupdate)
   {
       updateNum++;
       for(int i = 0; i < terminalList.size(); i ++)
       {
           for(int j = 0; j < sourceList.size(); j++)
           {
              negativeRes[i][j] = (-1)*flowTable[i][j];
              negaResDist[i][j] = (flowTable[i][j])? (-1)*distanceTable[i][j]: 0;
           }
       }

       for(int i = 0; i < terminalList.size(); i++)
       {
           for(int j = 0; j < sourceList.size(); j++)
           {
              positiveRes[i][j] = min((-1)*terminalList[i].maxflow, sourceList[j].maxflow)-flowTable[i][j];
              posResDist[i][j] = (positiveRes[i][j])? distanceTable[i][j]: 0;
           }
       }

       //Bellman's algorithm:
       //initialize:
       for(int j = 0; j < sourceList.size(); j++)
       {
           sourceList[j].minDist = 0;
           sourceList[j].NodeIndex = -1;
           sourceList[j].updated = false;
           sourceList[j].vistited = false;
           sourceList[j].NodeIndex = -1;
       }
       for(int i = 0; i < terminalList.size(); i++)
       {
           terminalList[i].minDist = infinte;
           terminalList[i].NodeIndex = -1;
           terminalList[i].updated = false;
           terminalList[i].vistited = false;
           terminalList[i].NodeIndex = -1;
       }
       trace.clear();
       //run algorithm
       bool flag = true;
       int iteration = 0;
       //while there being one node updated, continuing running the algorithm
       while(flag && (iteration < sourceList.size() + terminalList.size()))
       {     
          //initialize:
          for(int j = 0; j < sourceList.size(); j++) sourceList[j].updated = false;
          for(int i = 0; i < terminalList.size(); i++) terminalList[i].updated = false;
          flag = false;
          //start finding minumum distance
          for(int j = 0; j < sourceList.size(); j++)
          {
            for(int i = 0; i < terminalList.size(); i++)
            {
               if(terminalList[i].minDist == infinte) continue;
               else if(negaResDist[i][j])
               {
                  if(terminalList[i].minDist + negaResDist[i][j] < sourceList[j].minDist)
                  {
                    sourceList[j].minDist = terminalList[i].minDist + negaResDist[i][j];
                    sourceList[j].NodeIndex = i;
                    sourceList[j].updated = true;
                   }
                }
             }
           }
           for(int i = 0; i < terminalList.size(); i++)
           {
              for(int j = 0; j < sourceList.size(); j++)
              {
                if(posResDist[i][j])
                {
                  if(sourceList[j].minDist + posResDist[i][j] < terminalList[i].minDist)
                  {
                      terminalList[i].minDist = sourceList[j].minDist + posResDist[i][j];
                      terminalList[i].NodeIndex = j;
                      terminalList[i].updated = true;
                  }
                 }
               }
            }
            iteration++;
            for(int j = 0; j < sourceList.size(); j++)
            {
               flag = sourceList[j].updated;
               if(sourceList[j].updated) break;
            }
            for(int i = 0; i < terminalList.size(); i++)
            {
               flag = terminalList[i].updated;
               if(terminalList[i].updated) break;
            }
        }
        //(flag)? cout << "flag = true\n": cout << "flag = false\n";
        if(!flag) doupdate = false;
        else
        {
          //trace repeated nodes:
          bool repeated = false;
          for(int j = 0; j < sourceList.size(); j++)
          {
             if(sourceList[j].updated)
             {
                int curIndex = j;
                trace.push_back(curIndex);
                sourceList[j].vistited = true;
                curIndex = sourceList[j].NodeIndex;
                int counter = 1;
                while(!repeated && curIndex != -1)
                {
                   if(counter % 2)
                   {
                      repeated = terminalList[curIndex].vistited;
                      trace.push_back(curIndex);
                      terminalList[curIndex].vistited = true;
                      curIndex = terminalList[curIndex].NodeIndex; 
                   }
                   else if(counter % 2 == 0)
                   {
                      repeated = sourceList[curIndex].vistited;
                      trace.push_back(curIndex);
                      sourceList[curIndex].vistited = true;
                      curIndex = sourceList[curIndex].NodeIndex;            
                   }
                   counter++;
                }
                if(repeated) break;
             }
          }
          /*for(int i = 0; i < trace.size(); i++)
          {
            (i % 2)? cout << "t": cout << "s";
            cout << trace[i] << " ";
          }*/
          int traceSize = trace.size();
          int beginIndex = trace[traceSize - 1];
          int minUpdateFlow = infinte;
          int k = traceSize - 1;
          vector<int> cycleEdgeflow;
          cycleEdgeflow.clear();
          vector<int> cycleNode;
          cycleNode.clear();
          cycleNode.push_back(trace[k]);
          cycleNode.push_back(trace[k - 1]);
          if(traceSize % 2)
          {
            cycleEdgeflow.push_back(positiveRes[trace[k - 1]][trace[k]]);
            cycleEdgeflow.push_back(negativeRes[trace[k - 1]][trace[k - 2]]);
            cycleNode.push_back(trace[k]);
            cycleNode.push_back(trace[k - 1]);
          }
          else
          {
            cycleEdgeflow.push_back(negativeRes[trace[k]][trace[k - 1]]);
            cycleEdgeflow.push_back(positiveRes[trace[k - 1]][trace[k - 2]]);
            cycleNode.push_back(trace[k]);
            cycleNode.push_back(trace[k - 1]);
          }
          k = k - 2;
          int cycleEnd;
          if(repeated)
          {
             while(k && trace[k] != beginIndex)
             {
                cycleNode.push_back(trace[k]);
                cycleNode.push_back(trace[k - 1]);
                if(traceSize % 2)
                {
                  cycleEdgeflow.push_back(positiveRes[trace[k - 1]][trace[k]]);
                  cycleEdgeflow.push_back(negativeRes[trace[k - 1]][trace[k - 2]]);
                }
                else
                {
                  cycleEdgeflow.push_back(negativeRes[trace[k]][trace[k - 1]]);
                  cycleEdgeflow.push_back(positiveRes[trace[k - 2]][trace[k - 1]]);
                }
                k -= 2;
             }
             cycleEnd = k;
          }

          int minFlow = infinte;
          for(int i = 0; i < cycleEdgeflow.size(); i++)
          {
            if(minFlow > abs(cycleEdgeflow[i])) minFlow = abs(cycleEdgeflow[i]);
          }
            if(!minFlow) break;
            //cout << "minFlow = " << minFlow << endl;
            //update original graph
            for(int i = traceSize - 1; i > cycleEnd; i--)
            {
               if((i + 1) % 2) //trace[i] is source
               {
                  flowTable[trace[i - 1]][trace[i]] = flowTable[trace[i - 1]][trace[i]] + minFlow;
               }
               else //trace[i] is sink
               {
                  flowTable[trace[i]][trace[i - 1]] = flowTable[trace[i]][trace[i - 1]] - minFlow;
               }
            }
        }
     //(doupdate)? cout << "update again\n": cout << "optimized!\n";
       //cout << "updateNum = " << updateNum << endl;
     }
//end while
   //result calculating
   
   int wireArea = 0;
   vector<int> output;
   for(int i = 0; i < terminalList.size(); i++)
   {
      for(int j = 0; j < sourceList.size(); j++)
      {
      	if(flowTable[i][j])
      		{
      			output.push_back(sourceList[j].coordinate.first);
      			output.push_back(sourceList[j].coordinate.second);
      			output.push_back(terminalList[j].coordinate.first);
      			output.push_back(terminalList[j].coordinate.second);
      			output.push_back(flowTable[i][j]);
      		}
      	wireArea = wireArea + flowTable[i][j]*distanceTable[i][j];
      }
   }
   //cout << wireArea << endl;
   ofs << wireArea << endl;

   
   for(int i = 0; i < output.size(); i++)
   {
      ofs << output[i] << " ";
   	  if(i % 5 == 4) ofs << endl;
   }
   return 0;
}

bool myStr2Int(const string& str, int& num)
{
   num = 0;
   size_t i = 0;
   int sign = 1;
   if (str[0] == '-') { sign = -1; i = 1; }
   bool valid = false;
   for (; i < str.size(); ++i) {
      if (isdigit(str[i])) {
         num *= 10;
         num += int(str[i] - '0');
         valid = true;
      }
      else return false;
   }
   num *= sign;
   return valid;
}

int abs(int num)
{
    return (num > 0)? num: -num;
}

int min(int& a, int& b)
{
   return (a < b)? a: b;
}
