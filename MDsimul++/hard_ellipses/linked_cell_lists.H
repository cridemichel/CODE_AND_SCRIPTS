#include<vector>
#include "../mdlib/pvector.H"
using namespace std;
template <class ntype, class particles>
class linked_cell_lists
{
  vector<particles>* part;
  pvector<int,-1> inCell[3];
  pvector<int,-1> cellList;
  int parnum;
  pvector<ntype,2> L, L2;// box size
  int cellsx, cellsy, cellsz;
  void insert_in_new_cell(int i)
    {
      int n;
      n = (inCell[2][i] * cellsy + inCell[1][i])*cellsx + 
        inCell[0][i] + part->size();
      /* Inserimento di evIdA nella nuova cella (head) */
      cellList[i] = cellList[n];
      cellList[n] = i;

    }
  void remove_from_current_cell(int i)
    {
      int n;
      n = (inCell[2][i] * cellsy + inCell[1][i])*cellsx + inCell[0][i]
        + parnum;

      while (cellList[n] != i) 
        n = cellList[n];
      /* Eliminazione di evIdA dalla lista della cella n-esima */
      cellList[n] = cellList[i];
    }
public:
  linked_cell_lists(vector<particles>* P, pvector<ntype,2>& boxsize)
    {
      int k;
      part = P;
      L = boxsize;
      L2 = (1.0/2.0)*boxsize;
      parnum=part->size();
      for (k=0; k < 3; k++)
        inCell[k].allocate(part);
      cellList.allocate(parnum+cellsx*cellsy*cellsz);
    }
  void update_LL(int n)
    {
      int cox, coy, coz, cx, cy, cz;
      cox=inCell[0][n];
      coy=inCell[1][n];
      coz=inCell[2][n];
      cx =  ((*part)[n].r[0] + 0.5*L[0]) * cellsx / L[0];
      cy =  ((*part)[n].r[1] + 0.5*L[1]) * cellsy / L[1];
      cz =  ((*part)[n].r[2] + 0.5*L[2]) * cellsz / L[2];
      if (cx!=cox || cy!=coy || cz!=coz)
        {
          remove_from_current_cell(n);
          inCell[0][n] = cx;
          inCell[1][n] = cy;
          inCell[2][n] = cz;
          insert_in_new_cell(n);
        }
    }
  void build(void)
    {
      int j, n;

      for (j = 0; j < cellsx*cellsy*cellsz + parnum; j++)
        cellList[j] = -1;
      /* -1 vuol dire che non c'è nessuna particella nella cella j-esima */

      for (n = 0; n < parnum; n++)
        {
          inCell[0][n] =  ((*part)[n].r[0] + L2[0]) * cellsx / L[0];
          inCell[1][n] =  ((*part)[n].r[1] + L2[1]) * cellsy / L[1];
          inCell[2][n] =  ((*part)[n].r[2] + L2[2]) * cellsz / L[2];
          j = (inCell[2][n]*cellsy + inCell[1][n])*cellsx + 
            inCell[0][n] + parnum;
          cellList[n] = cellList[j];
          cellList[j] = n;
        }
    }
};
