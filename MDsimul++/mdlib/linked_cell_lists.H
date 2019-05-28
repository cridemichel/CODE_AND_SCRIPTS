#include<vector>
#include "../mdlib/pvector.H"
using namespace std;
// particles is the particle type, e.g. hardell<double>, hardcyl<double>, etc.
// rectangles are used as preliminary check of overlap (by using the separating axis theorem)
template <class particles, class rectangles, bool boxes=true>
class linked_cell_lists_2d
{
  using ntype = typename particles::numtype;
  vector<particles>* part;
  vector<rectangles> *rect;
  pvector<int,-1> inCell[2];
  pvector<int,-1> cellList, cellListBak;
  int parnum;
  pvector<ntype,2> L, L2;// box size
  pvector<ntype,2> Lbak, L2bak;
  int cellsx, cellsy;
  int cellsxBak, cellsyBak;
  int parnumBak;
  pvector<int,-1> inCellBak[2];
  int cxini, cyini;
  ntype rcut;
  void insert_in_new_cell(int i)
    {
      int n;
      n = inCell[1][i]*cellsx + 
        inCell[0][i] + parnum;
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
      if constexpr (boxes==false)
        {
          if (overlap(A,B) < 0.0)
            {
              return true;
            }
        }
      else
        {
          if (overlap(rectA,rectB) < 0.0)
            {
              if (overlap(A,B) < 0.0)
                {
                  //cout << "true\n";
#if 0
                  if (na==40||n==40)
                    {
                      cout << "na=" << na << " n" << n << "\n";
                      A.show(">>>A");
                      B.show(">>>B");
                    }
#endif
                  return true;
                }
            }
        }
      //cout << "false\n";
      return false;
    }
public:
  linked_cell_lists_2d(vector<particles>* P, vector<rectangles> *R, pvector<ntype,2>& boxsize)
    {
      init(P,R,boxsize);    
    }
  linked_cell_lists_2d()
    {
      L[0]=L[1]=-1.0;
    }
  void init(vector<particles>* P, vector<rectangles> *R, pvector<ntype,2>& boxsize, ntype rc)
    {
      int k;
      part = P;
      rect = R;
      rcut = rc;
      L = boxsize;
      cellsx = L[0] / rcut;
      cellsy = L[1] / rcut;
      //cout << "L=" << L[0] << L[1] << " rcut=" << rcut << "\n";
      L2 = (1.0/2.0)*boxsize;
      parnum=part->size();
      //cout << "cellsx= " << cellsx << " cellsy= " << cellsy << "\n";
      for (k=0; k < 2; k++)
        {
          inCell[k].allocate(parnum);
          inCellBak[k].allocate(parnum);
        }
      cellList.allocate(parnum+cellsx*cellsy);
      cellListBak.allocate(parnum+cellsx*cellsy);
    }
  void set_ini_numcells(void)
    {
      cxini = cellsx;
      cyini = cellsy;
    }
  void update_numcells(pvector<ntype,2> Lb)
    {
      L = Lb;
      L2 = 0.5*L;
      cellsx = L[0] / rcut;
      cellsy = L[1] / rcut;
      if (cellsx*cellsy > cxini*cyini)
        {
          //cellList = realloc(cellList, sizeof(int)*(cellsx*cellsy*cellsz+Oparams.parnum));
          cellList.resize(parnum+cellsx*cellsy);
          cellListBak.resize(parnum+cellsx*cellsy);
        } 
    }
  void store(void)
    {
      parnumBak=parnum;
      Lbak=L;
      L2bak=L2;
      cellsxBak = cellsx;
      cellsyBak = cellsy;
      for (auto k=0; k < parnum + cellsx*cellsy; k++)
        {
          cellListBak[k] = cellList[k];
        }
      for (auto i=0; i < parnum; i++)
        {
          inCellBak[0][i] = inCell[0][i];
          inCellBak[1][i] = inCell[1][i];
        }
    }
  void restore(void)
    {
      L2=L2bak;
      L=Lbak;
      cellsx = cellsxBak;
      cellsy = cellsyBak; 
      for (auto k=0; k < parnum + cellsx*cellsy; k++)
        {
          cellList[k] = cellListBak[k];
        }
      for (auto i=0; i < parnum; i++)
        {
          inCell[0][i] = inCellBak[0][i];
          inCell[1][i] = inCellBak[1][i];
        }
    }
  bool check_overlap(int ip, int nb=-1)
    {
      int kk, iX, jX, iY, jY, n, na;
      int cellRange[4];
      pvector<ntype,2> shift;
      na=ip;
      for (kk = 0;  kk < 2; kk++)
        {
          cellRange[2*kk]   = - 1;
          cellRange[2*kk+1] =   1;
        }

      for (iY = cellRange[2]; iY <= cellRange[3]; iY++) 
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
          for (iX = cellRange[0]; iX <= cellRange[1]; iX++) 
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
                          //cout << "overlap found na=" << na << " n=" << n << "\n";
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
