#include "iostream"
using namespace std;

// Adjust Q and R to get the more accurate reading
// Q : Kalman Filter model process noise
// R: Sensor noise , adjust this value to minimize the effect of noise in sensor's reading
class KalmanFilter {
  public:
    int n; // State dimension
    int m; // Control dimension
    int k; // Measurement dimension

    float* s; // State vector [n]
    float* P; // Covariance matrix [n*n]
    float* A; // Transition matrix [n*n]
    float* B; // Control matrix [n*m]
    float* H; // Measurement matrix [k*n]
    float* Q; // Process noise covariance [n*n]
    float* R; // Measurement noise covariance [k*k]
    float* I; // Identity matrix [n*n]

    KalmanFilter(int state_dim, int control_dim, int measure_dim) {
      n = state_dim; m = control_dim; k = measure_dim;
      s = new float[n];
      P = new float[n * n];
      A = new float[n * n];
      B = new float[n * m];
      H = new float[k * n];
      Q = new float[n * n];
      R = new float[k * k];
      I = new float[n * n];
      
      // Initialize Identity Matrix
      for(int i=0; i<n*n; i++) I[i] = (i % (n+1) == 0) ? 1.0 : 0.0;
    }

    void predict(float* u) {
      // 1. st = A*st-1 + B*ut
      float* new_s = new float[n];
      matrixMultiply(A, s, n, n, 1, new_s);
      if (u != nullptr) {
        float* Bu = new float[n];
        matrixMultiply(B, u, n, m, 1, Bu);
        for(int i=0; i<n; i++) new_s[i] += Bu[i];
        delete[] Bu;
      }
      for(int i=0; i<n; i++) s[i] = new_s[i];
      delete[] new_s;

      // 2. Pt = A*Pt-1*A^T + Q
      float* AT = new float[n * n];
      float* AP = new float[n * n];
      float* APAT = new float[n * n];
      transpose(A, AT, n, n);
      matrixMultiply(A, P, n, n, n, AP);
      matrixMultiply(AP, AT, n, n, n, APAT);
      for(int i=0; i<n*n; i++) P[i] = APAT[i] + Q[i];
      
      delete[] AT; delete[] AP; delete[] APAT;
    }

    void update(float* z) {
      // 3. K = PH^T * (HPH^T + R)^-1
      float* HT = new float[n * k];
      float* PH_T = new float[n * k];
      float* HPH_T = new float[k * k];
      float* S = new float[k * k]; // Innovation covariance
      transpose(H, HT, k, n);
      matrixMultiply(P, HT, n, n, k, PH_T);
      matrixMultiply(H, PH_T, k, n, k, HPH_T);
      for(int i=0; i<k*k; i++) S[i] = HPH_T[i] + R[i];

      float* S_inv = new float[k * k];
      invertMatrix(S, S_inv, k);

      float* K = new float[n * k];
      matrixMultiply(PH_T, S_inv, n, k, k, K);

      // 4. st = st + K(zt - Hst)
      float* Hs = new float[k];
      matrixMultiply(H, s, k, n, 1, Hs);
      float* innovation = new float[k];
      for(int i=0; i<k; i++) innovation[i] = z[i] - Hs[i];
      float* K_innov = new float[n];
      matrixMultiply(K, innovation, n, k, 1, K_innov);
      for(int i=0; i<n; i++) s[i] += K_innov[i];

      // 5. Pt = (I - KH)P
      float* KH = new float[n * n];
      float* I_KH = new float[n * n];
      matrixMultiply(K, H, n, k, n, KH);
      for(int i=0; i<n*n; i++) I_KH[i] = I[i] - KH[i];
      float* new_P = new float[n * n];
      matrixMultiply(I_KH, P, n, n, n, new_P);
      for(int i=0; i<n*n; i++) P[i] = new_P[i];

      // Cleanup
      delete[] HT; delete[] PH_T; delete[] HPH_T; delete[] S; 
      delete[] S_inv; delete[] K; delete[] Hs; delete[] innovation;
      delete[] K_innov; delete[] KH; delete[] I_KH; delete[] new_P;
    }

    void init(float* target, float* source, int r, int c) {
        for (int i = 0; i < r * c; i++) {
            target[i] = source[i];
        }
    }

  private:
    void matrixMultiply(float* A, float* B, int r1, int c1, int c2, float* res) {
      for (int i = 0; i < r1; i++) {
        for (int j = 0; j < c2; j++) {
          res[i * c2 + j] = 0;
          for (int l = 0; l < c1; l++) res[i * c2 + j] += A[i * c1 + l] * B[l * c2 + j];
        }
      }
    }

    void transpose(float* A, float* res, int r, int c) {
      for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) res[j * r + i] = A[i * c + j];
      }
    }

  void invertMatrix(float* A, float* res, int n) {
    if (n == 1) {
      res[0] = 1.0 / A[0];
    } else if (n == 2) {
      float det = A[0] * A[3] - A[1] * A[2];
      if (det == 0) return; // Add check to prevent division by zero
      res[0] = A[3] / det; res[1] = -A[1] / det;
      res[2] = -A[2] / det; res[3] = A[0] / det;
    } else if (n == 3) {
      // 3x3 Matrix Inversion for Position, Velocity, Acceleration states
      // Matrix indices: 
      // 0 1 2
      // 3 4 5
      // 6 7 8
      
      float det = A[0] * (A[4] * A[8] - A[5] * A[7]) -
                  A[1] * (A[3] * A[8] - A[5] * A[6]) +
                  A[2] * (A[3] * A[7] - A[4] * A[6]);

      if (det == 0) return; // Ensure the matrix is non-singular
      float invDet = 1.0 / det;

      res[0] = (A[4] * A[8] - A[5] * A[7]) * invDet;
      res[1] = (A[2] * A[7] - A[1] * A[8]) * invDet;
      res[2] = (A[1] * A[5] - A[2] * A[4]) * invDet;
      res[3] = (A[5] * A[6] - A[3] * A[8]) * invDet;
      res[4] = (A[0] * A[8] - A[2] * A[6]) * invDet;
      res[5] = (A[2] * A[3] - A[0] * A[5]) * invDet;
      res[6] = (A[3] * A[7] - A[4] * A[6]) * invDet;
      res[7] = (A[1] * A[6] - A[0] * A[7]) * invDet;
      res[8] = (A[0] * A[4] - A[1] * A[3]) * invDet;
    }
  }
};

int main() {
    float dt = 0.1; // Sampling time as used in the paper
    
    // Initialize n=3 (d,v,a), m=0 (no control), k=1 (acceleration measure)
    KalmanFilter kf(3, 0, 1);

    // 1. Transition Matrix A (3x3)
    // d = d + v*dt + 0.5*a*dt^2
    // v = v + a*dt
    // a = a
    float a_data[] = {
        1,  dt,  0.5f*dt*dt,
        0,  1,   dt,
        0,  0,   1
    };
    kf.init(kf.A, a_data, 3, 3);

    // 2. Measurement Matrix H (1x3)
    // We only measure 'a' (the 3rd element of the state)
    float h_data[] = {0, 0, 1};
    kf.init(kf.H, h_data, 1, 3);

    // 3. Process Noise Q (3x3) - Small values for stability
    float q_data[] = {
        0.01, 0,    0,
        0,    0.01, 0,
        0,    0,    0.01
    };
    kf.init(kf.Q, q_data, 3, 3);

    // 4. Measurement Noise R (1x1) - Variance of your accelerometer
    float r_data[] = {0.1}; 
    kf.init(kf.R, r_data, 1, 1);

    // 5. Initial Covariance P (3x3) - Identity or high uncertainty
    float p_data[] = {
        1, 0, 0,
        0, 1, 0,
        0, 0, 1
    };
    kf.init(kf.P, p_data, 3, 3);

    float test_data[] = {2,1,1.5,3}; //accelerometer data of one axis
    int count=0;
    float s_o[] = {0, 12, 1};
    kf.init(kf.s, s_o, 3, 1);
    while(count<4){
      kf.predict(nullptr);
      cout<<"phase: "<<count<<endl;
      cout<<"Predicted: "<<endl;
      for(int i=0;i<3;i++){
        cout<< kf.s[i]<<endl;
      }
      kf.update(&test_data[count]);
      cout<<"Corrected: "<<endl;
      for(int i=0;i<3;i++){
        cout<< kf.s[i]<<endl;
      }
      cout<<endl;
      count++;
    }
    return 0;
}