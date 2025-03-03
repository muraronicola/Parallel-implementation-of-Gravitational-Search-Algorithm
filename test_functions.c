
float sphere(float* x, int dim) {
    float sum = 0;
    for (int i = 0; i < dim; i++) {
        sum += x[i] * x[i];
    }
    return sum;
}