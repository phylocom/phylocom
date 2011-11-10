
void meanvar(float X[], unsigned int n, float *mean, float *var);
void ttest(float X[], unsigned int nx, float Y[], unsigned int ny, \
           float *t, float *p);
void correl(float X[], float Y[], unsigned int n, \
            float *r, float *p);
void contabg(int **X, unsigned int nx, \
             unsigned int ny, float *g, float *p);
float ttable(float score, int dfin, int tailed);
float chitable(float score, int dfin);
