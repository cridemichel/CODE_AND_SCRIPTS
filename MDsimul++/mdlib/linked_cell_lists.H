#include<vector>
#include "../mdlib/pvector.H"
using namespace std;
// particles is the particle type, e.g. hardell<double>, hardcyl<double>, etc.
// rectangles are used as preliminary check of overlap (by using the separating axis theorem)
template <class particles, class rectangles>
class linked_cell_lists_2d
{
  using ntype = typename particles::numtype;
  vector<particles>* part;
  vector<rectangles> *rect;
  pvector<int,-1> inCell[2];
  pvector<int,-1> cellList;
  int parnum;
  pvector<ntype,2> L, L2;// box size
  int cellsx, cellsy;
  void insert_in_new_cell(int i)
    {
      int n;
      n = inCell[1][i]*cellsx + 
        inCell[0][i] + part->size();
      /* Inserimento di evIdA nella nuova cella (head) */
      cellList[i] = cellList[n];
      cellList[n] = i;

    }
  void remove_from_current_cell(int i)
    {
      int n;
      n = inCell[1][i]*cellsx + inCell[0][i]
        + parnum;

      while (cellList[n] != i) 
        n = cellList[n];
      /* Eliminazione di evIdA dalla lista della cella n-esima */
      cellList[n] = cellList[i];
    }
  
  bool overlap_two(int na, int n, pvector<ntype,2> shift)
    {
      auto rectA = (*rect)[na];
      auto rectB = (*rect)[n];
      auto A = (*part)[na];
      auto B = (*part)[n];
      rectA.shift << 0,0;
      rectB.shift = shift;
      A.shift << 0, 0;
      B.shift = shift; 
      if (overlap(rectA,rectB) < 0.0)
        {
          if (overlap(A,B) < 0.0)
            return true;
        }
      return false;
    }
public:
  linked_cell_lists_2d(vector<particles>* P, pvector<ntype,2>& boxsize)
    {
      init(P,boxsize);    
    }
  linked_cell_lists_2d()
    {
      L[0]=L[1]=-1.0;
    }
  void init(vector<particles>* P, pvector<ntype,2>& boxsize)
    {
      int k;
      part = P;
      L = boxsize;
      L2 = (1.0/2.0)*boxsize;
      parnum=part->size();
      for (k=0; k < 3; k++)
        inCell[k].allocate(parnum);
      cellList.allocate(parnum+cellsx*cellsy);
    }
  
  bool check_overlap(int ip)
    {
      int kk, nb, iX, jX, iY, jY, n, na;
      int cellRange[4];
      pvector<ntype,2> shift;
      na=ip;
      nb=-1;
      for (kk = 0;  kk < 2; kk++)
        {
          cellRange[2*kk]   = - 1;
          cellRange[2*kk+1] =   1;
        }

      for (iY = cellRange[2]; iY <= cellRange[3]; iY ++) 
        {
          jY = inCell[1][na] + iY;    
          shift[1] = 0.0;
          if (jY == -1) 
            {
              jY = cellsy - 1;    
              shift[1] = -L[1];
            } 
          else if (jY == cellsy) 
            {
              jY = 0;    
              shift[1] = L[1];
            }
          for (iX = cellRange[0]; iX <= cellRange[1]; iX ++) 
            {
              jX = inCell[0][na] + iX;    
              shift[0] = 0.0;
              if (jX == -1) 
                {
                  jX = cellsx - 1;    
                  shift[0] = - L[0];
                } 
              else if (jX == cellsx) 
                {
                  jX = 0;   
                  shift[0] = L[0];
                }
              n = jY*cellsx + jX + parnum;
              for (n = cellList[n]; n > -1; n = cellList[n]) 
                {
                  if (n != na && n != nb && (nb >= -1 || n < na)) 
                    {
                      if (overlap_two(na, n, shift)==true)
                        {
                          return true;
                        }
                    }
                } 
            }
        }
      return false;
    }
  void update_LL(int n)
    {
      int cox, coy, cx, cy;
      cox=inCell[0][n];
      coy=inCell[1][n];
      cx =  ((*part)[n].r[0] + 0.5*L[0]) * cellsx / L[0];
      cy =  ((*part)[n].r[1] + 0.5*L[1]) * cellsy / L[1];
      if (cx!=cox || cy!=coy)
        {
          remove_from_current_cell(n);
          inCell[0][n] = cx;
          inCell[1][n] = cy;
          insert_in_new_cell(n);
        }
    }
  void build(void)
    {
      int j, n;

      for (j = 0; j < cellsx*cellsy + parnum; j++)
        cellList[j] = -1;
      /* -1 vuol dire che non c'è nessuna particella nella cella j-esima */

      for (n = 0; n < parnum; n++)
        {
          inCell[0][n] =  ((*part)[n].r[0] + L2[0]) * cellsx / L[0];
          inCell[1][n] =  ((*part)[n].r[1] + L2[1]) * cellsy / L[1];
          j = inCell[1][n]*cellsx + 
            inCell[0][n] + parnum;
          cellList[n] = cellList[j];
          cellList[j] = n;
        }
    }
};
