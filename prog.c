
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

// вычисление скалярного произведения
double scalar_product(double **u, double **v, double h1, double h2, int M,
                      int N) {
  double ans = 0.0;
  for (int j = 1; j < N - 1; ++j) {
    for (int i = 1; i < M - 1; ++i) {
      ans += h1 * h2 * u[j][i] * v[j][i];
    }
  }
  return ans;
}

// вычисление нормы L2
double norm(double **u, double h1, double h2, int M, int N) {
  return sqrt(scalar_product(u, u, h1, h2, M, N));
}

// копирование матрицы
void mat_copy(double **src, double **target, int M, int N) {
  for (int j = 0; j < N; ++j) {
    for (int i = 0; i < M; ++i) {
      target[j][i] = src[j][i];
    }
  }
}

// установка во все значения матрицы значения val
void mat_set_value(double **u, int M, int N, int val) {
  for (int j = 0; j < N; ++j) {
    for (int i = 0; i < M; ++i) {
      u[j][i] = val;
    }
  }
}

// создние матрицы и инициализация всех ее элементов нулем
double **mat_create(int M, int N) {
  double **mat = (double **)calloc(N, sizeof(double *));
  for (int j = 0; j < N; ++j) {
    mat[j] = (double *)calloc(M, sizeof(double));
  }
  mat_set_value(mat, M, N, 0);
  return mat;
}

// освобождение памяти выделенной под матрицу
void mat_free(double **mat, int M, int N) {
  for (int j = 0; j < N; ++j) {
    free(mat[j]);
  }
  free(mat);
}

// максимум двую чисел
double max(double a, double b) {
  if (a > b)
    return a;
  return b;
}

// вычитание из матицы u матрицы v поэлементно, создаем новую матрицу !!!
void mat_minus(double **u, double **v, int M, int N, double **ans) {
  for (int j = 1; j < N - 1; ++j) {
    for (int i = 1; i < M - 1; ++i) {
      ans[j][i] = u[j][i] - v[j][i];
    }
  }
  return;
}

// оператор A
void A_fun(double **a, double **b, double **w, int M, int N, double h1,
           double h2, double **ans) {
  for (int j = 1; j < N; ++j) {
    for (int i = 1; i < M; ++i) {
      double s1 = ((a[j][i + 1] * (w[j][i + 1] - w[j][i]) / h1) -
                   (a[j][i] * (w[j][i] - w[j][i - 1]) / h1)) /
                  h1;
      double s2 = ((b[j + 1][i] * (w[j + 1][i] - w[j][i]) / h2) -
                   (b[j][i] * (w[j][i] - w[j - 1][i]) / h2)) /
                  h2;
      ans[j][i] = -s1 - s2;
    }
  }
  return;
}

// оператор B
void B_fun(double **F, int M, int N, double **ans) {
  for (int j = 1; j < N; ++j) {
    for (int i = 1; i < M; ++i) {
      ans[j][i] = F[j][i];
    }
  }
  return;
}

// умножение всех элементов матрицы на число val
void mat_mul_number(double **mat, double val, int M, int N) {
  for (int j = 1; j < N; ++j) {
    for (int i = 1; i < M; ++i) {
      mat[j][i] *= val;
    }
  }
  return;
}

// структура реализующая точку плоскости
struct Point {
  double x;
  double y;
} typedef Point;

// вычисление l_ij
double calc_l_ij(Point P1, Point P2) {
  if (P1.x != P2.x) {
    fprintf(stderr, "Error in calc_l_ij, P1.x != P2.x\n");
    return 0;
  }
  if (P1.y >= P2.y) {
    fprintf(stderr, "Error in calc_l_ij, P1.y >= P2.y\n");
    return 0;
  }
  double x = P1.x;
  double y = (-4.0 / 3.0) * x + 4;
  if (y > P2.y) {
    return P2.y - P1.y;
  } else if (y < P1.y) {
    return 0;
  } else {
    return y - P1.y;
  }
  fprintf(stderr, "Error in calc_l_ij, default return\n");
  return 0;
}

// вычисление a_ij
double calc_a_ij(double h2, Point p1, Point p2, double eps) {
  double l_ij = calc_l_ij(p1, p2);
  double a_ij = (l_ij / h2) + (1 - l_ij / h2) / eps;
  return a_ij;
}

void test_calc_l_ij() {
  {
    Point p1;
    p1.x = 1;
    p1.y = 1;
    Point p2;
    p2.x = 1;
    p2.y = 2;
    double ans = calc_l_ij(p1, p2);
    if (-0.0001 < (ans - 1.0) && (ans - 1.0) < 0.0001) {
      printf("Pass test1 for calc_l_ij\n");
    } else {
      printf("Faild test1 for calc_l_ij\n");
    }
  }
  {
    Point p1;
    p1.x = 1;
    p1.y = 2;
    Point p2;
    p2.x = 1;
    p2.y = 3;
    double ans = calc_l_ij(p1, p2);
    if (-0.0001 < (ans - (2.0 / 3.0)) && (ans - (2.0 / 3.0)) < 0.0001) {
      printf("Pass test2 for calc_l_ij\n");
    } else {
      printf("Faild test2 for calc_l_ij\n");
    }
  }
  {
    Point p1;
    p1.x = 1;
    p1.y = 3;
    Point p2;
    p2.x = 1;
    p2.y = 4;
    double ans = calc_l_ij(p1, p2);
    if (-0.0001 < (ans - 0.0) && (ans - 0.0) < 0.0001) {
      printf("Pass test3 for calc_l_ij\n");
    } else {
      printf("Faild test3 for calc_l_ij\n");
    }
  }
}

// вычисление p_ij
double calc_p_ij(Point P1, Point P2) {
  if (P1.y != P2.y) {
    fprintf(stderr, "Error in calc_p_ij, P1.y != P2.y\n");
    return 0;
  }
  if (P1.x >= P2.x) {
    fprintf(stderr, "Error in calc_p_ij, P1.x >= P2.x\n");
    return 0;
  }
  double y = P1.y;
  double x = (-3.0 / 4.0) * y + 3;
  if (x > P2.x) {
    return P2.x - P1.x;
  } else if (x < P1.x) {
    return 0;
  } else {
    return x - P1.x;
  }
  fprintf(stderr, "Error in calc_p_ij, default return\n");
  return 0;
}

// вычисление b_ij
double calc_b_ij(double h1, Point p1, Point p2, double eps) {
  double p_ij = calc_p_ij(p1, p2);
  double b_ij = (p_ij / h1) + (1 - p_ij / h1) / eps;
  return b_ij;
}

void test_calc_p_ij() {
  {
    Point p1;
    p1.x = 1;
    p1.y = 1;
    Point p2;
    p2.x = 2;
    p2.y = 1;
    double ans = calc_p_ij(p1, p2);
    if (-0.0001 < (ans - 1.0) && (ans - 1.0) < 0.0001) {
      printf("Pass test1 for calc_p_ij\n");
    } else {
      printf("Faild test1 for calc_p_ij\n");
    }
  }
  {
    Point p1;
    p1.x = 2;
    p1.y = 1;
    Point p2;
    p2.x = 3;
    p2.y = 1;
    double ans = calc_p_ij(p1, p2);
    if (-0.0001 < (ans - (1.0 / 4.0)) && (ans - (1.0 / 4.0)) < 0.0001) {
      printf("Pass test2 for calc_p_ij\n");
    } else {
      printf("Faild test2 for calc_p_ij\n");
    }
  }
  {
    Point p1;
    p1.x = 3;
    p1.y = 1;
    Point p2;
    p2.x = 4;
    p2.y = 1;
    double ans = calc_p_ij(p1, p2);
    if (-0.0001 < (ans - 0.0) && (ans - 0.0) < 0.0001) {
      printf("Pass test3 for calc_p_ij\n");
    } else {
      printf("Faild test3 for calc_p_ij\n");
    }
  }
  {
    Point p1;
    p1.x = 1;
    p1.y = 2;
    Point p2;
    p2.x = 2;
    p2.y = 2;
    double ans = calc_p_ij(p1, p2);
    if (-0.0001 < (ans - 0.5) && (ans - 0.5) < 0.0001) {
      printf("Pass test4 for calc_p_ij\n");
    } else {
      printf("Faild test4 for calc_p_ij\n");
    }
  }
}

// вычисление S_ij
double calc_S_ij(double h1, double h2, Point p) {
  Point p1, p2, p3, p4;
  p1.x = p.x - h1 / 2;
  p1.y = p.y - h2 / 2;
  p2.x = p.x - h1 / 2;
  p2.y = p.y + h2 / 2;
  p3.x = p.x + h1 / 2;
  p3.y = p.y + h2 / 2;
  p4.x = p.x + h1 / 2;
  p4.y = p.y - h2 / 2;
  {
    // прямоугольник лежит выше
    double x = p1.x;
    double y = (-4.0 / 3.0) * x + 4.0;
    if (y <= p1.y) {
      return 0;
    }
  }
  {
    // прямоугольник лежит ниже
    double x = p3.x;
    double y = (-4.0 / 3.0) * x + 4.0;
    if (y >= p3.y) {
      return h1 * h2;
    }
  }
  {
    Point p_l, p_r, p_u, p_d;
    p_l.x = p1.x;
    p_l.y = -(4.0 / 3.0) * p_l.x + 4.0;
    p_r.x = p3.x;
    p_r.y = -(4.0 / 3.0) * p_r.x + 4.0;
    p_d.y = p1.y;
    p_d.x = -(3.0 / 4.0) * p_d.y + 3.0;
    p_u.y = p3.y;
    p_u.x = -(3.0 / 4.0) * p_u.y + 3.0;

    if (p1.y <= p_l.y && p_l.y <= p2.y) {
      if (p1.x <= p_d.x && p_d.x <= p4.x) {
        return 0.5 * (p_l.y - p1.y) * (p_d.x - p1.x);
      } else {
        return 0.5 * (p4.x - p1.x) * ((p_l.y - p1.y) + (p_r.y - p4.y));
      }
    } else {
      if (p1.x <= p_d.x && p_d.x <= p4.x) {
        return 0.5 * (p2.y - p1.y) * ((p_u.x - p2.x) + (p_d.x - p1.x));
      } else {
        return h1 * h2 - 0.5 * (p3.x - p_u.x) * (p3.y - p_r.y);
      }
    }
  }
  fprintf(stderr, "Error in calc_p_ij, default return\n");
  return 0;
}

void test_calc_S_ij() {
  {
    double h1 = 0.3, h2 = 0.4;
    Point p;
    p.x = 1.0, p.y = 1.0;
    double ans = calc_S_ij(h1, h2, p);
    if (-0.0001 < (ans - (0.3 * 0.4)) && (ans - (0.3 * 0.4)) < 0.0001) {
      printf("Pass test1 for calc_S_ij\n");
    } else {
      printf("Faild test1 for calc_S_ij\n");
    }
  }
  {
    double h1 = 0.3, h2 = 0.4;
    Point p;
    p.x = 1.0, p.y = 2.5;
    double ans = calc_S_ij(h1, h2, p);
    if (-0.0001 < (ans - (0.0995)) && (ans - (0.0995)) < 0.0001) {
      printf("Pass test2 for calc_S_ij\n");
    } else {
      printf("Faild test2 for calc_S_ij\n");
    }
  }
  {
    double h1 = 0.3, h2 = 0.4;
    Point p;
    p.x = 1.5, p.y = 2.0;
    double ans = calc_S_ij(h1, h2, p);
    if (-0.0001 < (ans - (0.06)) && (ans - (0.06)) < 0.0001) {
      printf("Pass test3 for calc_S_ij\n");
    } else {
      printf("Faild test3 for calc_S_ij\n");
    }
  }
  {
    double h1 = 0.3, h2 = 0.4;
    Point p;
    p.x = 2.5, p.y = 1.0;
    double ans = calc_S_ij(h1, h2, p);
    if (-0.0001 < (ans - (0.0016)) && (ans - (0.0016)) < 0.0001) {
      printf("Pass test4 for calc_S_ij\n");
    } else {
      printf("Faild test4 for calc_S_ij\n");
    }
  }
  {
    double h1 = 0.3, h2 = 0.4;
    Point p;
    p.x = 2.5, p.y = 3.0;
    double ans = calc_S_ij(h1, h2, p);
    if (-0.0001 < (ans - (0.0)) && (ans - (0.0)) < 0.0001) {
      printf("Pass test5 for calc_S_ij\n");
    } else {
      printf("Faild test5 for calc_S_ij\n");
    }
  }
}

// вычисление F_ij
double calc_F_ij(double h1, double h2, Point p) {
  double S = calc_S_ij(h1, h2, p);
  double F_ij = S / (h1 * h2); // f(x_i, y_i) == 1;
  return F_ij;
}

// инициализируем и заполняем матрицу a
double **init_a(int M, int N, double h1, double h2, double eps) {
  double **a = mat_create(M + 1, N + 1);
  for (int j = 1; j < N + 1; ++j) {
    for (int i = 1; i < M + 1; ++i) {
      Point p1;
      p1.x = i * h1 - 0.5 * h1, p1.y = j * h2 - 0.5 * h2;
      Point p2;
      p2.x = i * h1 - 0.5 * h1, p2.y = j * h2 + 0.5 * h2;
      a[j][i] = calc_a_ij(h2, p1, p2, eps);
    }
  }
  return a;
}

// инициализируем и заполняем матрицу b
double **init_b(int M, int N, double h1, double h2, double eps) {
  double **b = mat_create(M + 1, N + 1);
  for (int j = 1; j < N + 1; ++j) {
    for (int i = 1; i < M + 1; ++i) {
      Point p1;
      p1.x = i * h1 - 0.5 * h1, p1.y = j * h2 - 0.5 * h2;
      Point p2;
      p2.x = i * h1 + 0.5 * h1, p2.y = j * h2 - 0.5 * h2;
      b[j][i] = calc_b_ij(h1, p1, p2, eps);
    }
  }
  return b;
}

// инициализируем и заполняем матрицу F
double **init_F(int M, int N, double h1, double h2, double eps) {
  double **F = mat_create(M, N);
  for (int j = 1; j < N; ++j) {
    for (int i = 1; i < M; ++i) {
      Point p;
      p.x = i * h1, p.y = j * h2;
      F[j][i] = calc_F_ij(h1, h2, p);
    }
  }
  return F;
}

// вывд матрицы в точностью 5 знаков после запятой
void mat_print(double **mat, int M, int N) {
  for (int j = N - 1; j >= 0; --j) {
    for (int i = 0; i < M; ++i) {
      printf("%.5f, ", mat[j][i]);
    }
    printf("\n");
  }
  printf("\n\n\n");
}

int main() {
  // запуск тестов
  if (0) {
    test_calc_l_ij();
    test_calc_p_ij();
    test_calc_S_ij();
  }
  double start, end;
  start = omp_get_wtime();

  const double A1 = 0.0, A2 = 0.0, B1 = 3.0, B2 = 4.0;
  const int N = 40, M = 40;
  const double h1 = (B1 - A1) / (1.0 * M);
  const double h2 = (B2 - A2) / (1.0 * N);
  const int count_iter = 10;
  const double delta = 0.000001;
  const double eps = max(h1, h2) * max(h1, h2);

  double **a, **b, **F;
  a = init_a(M, N, h1, h2, eps);
  b = init_b(M, N, h1, h2, eps);
  F = init_F(M, N, h1, h2, eps);

  double **w_k = mat_create(M + 1, N + 1);
  double **w_k_plus1 = mat_create(M + 1, N + 1);
  double **temp1 = mat_create(M, N);
  double **temp2 = mat_create(M, N);
  double **temp3 = mat_create(M, N);
  double **r_k = mat_create(M + 1, N + 1);

  printf("Start iteration\n\n");

  int i = 0;
  double err = delta + 1;
  for (; err > delta; ++i) {
    A_fun(a, b, w_k, M, N, h1, h2, temp1);
    B_fun(F, M, N, temp2);
    mat_minus(temp1, temp2, M + 1, N + 1, r_k);

    A_fun(a, b, r_k, M, N, h1, h2, temp3);
    double tau_k_plus1 = scalar_product(r_k, r_k, h1, h2, M, N) /
                         scalar_product(temp3, r_k, h1, h2, M, N);

    mat_mul_number(r_k, tau_k_plus1, M, N);

    err = norm(r_k, h1, h2, M, N);
    // printf("%f, ", err);
    mat_minus(w_k, r_k, M + 1, N + 1, w_k_plus1);

    mat_copy(w_k_plus1, w_k, M + 1, N + 1);
    // mat_print(w_k, M+1, N+1);
  }

  printf("\n");
  printf("Count iteration: %i\n", i);
  // printf("Result w:\n");
  // mat_print(w_k, M+1, N+1);

  mat_free(a, M + 1, N + 1);
  mat_free(b, M + 1, N + 1);
  mat_free(F, M, N);
  mat_free(temp1, M, N);
  mat_free(temp2, M, N);
  mat_free(temp3, M, N);
  mat_free(w_k, M + 1, N + 1);
  mat_free(w_k_plus1, M + 1, N + 1);
  mat_free(r_k, M + 1, N + 1);

  end = omp_get_wtime();
  printf("Work took %f seconds\n", end - start);
  return 0;
}
