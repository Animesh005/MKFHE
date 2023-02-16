#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

const int INPUT_WIDTH = 28;
const int INPUT_HEIGHT = 28;
const int NUM_CLASSES = 10;
const int WINDOW_SIZE = 4;
const int FILTER_SIZE = 4;
const int NUM_FILTERS = 5;
const int STRIDE = 2;
const int NUM_TRAIN_IMAGES = 60000;
const int NUM_TRAIN_LABELS = 60000;
const int NUM_TEST_IMAGES = 10000;
const int NUM_TEST_LABELS = 10000;

vector<vector<vector<int>>> load_train_data()
{
    std::ifstream file_train_data("test/mnist_data/train-images.idx3-ubyte", std::ios::binary);

    if (!file_train_data.is_open()) {
        std::cerr << "Failed to open image file\n";
    }

    // Skip the first 16 bytes of the file, which contain header information
    file_train_data.seekg(16);

    // Read the contents of the file into a vector
    std::vector<unsigned char> train_data((std::istreambuf_iterator<char>(file_train_data)),
                                    std::istreambuf_iterator<char>());

    // Create a 3D array to store the image data
    std::vector<std::vector<std::vector<int>>> images(NUM_TRAIN_IMAGES, std::vector<std::vector<int>>(INPUT_HEIGHT, std::vector<int>(INPUT_WIDTH)));

    // Fill the image data from the vector
    for (int i = 0; i < NUM_TRAIN_IMAGES; i++) {
        for (int j = 0; j < INPUT_HEIGHT; j++) {
            for (int k = 0; k < INPUT_WIDTH; k++) {
                images[i][j][k] = ((int)train_data[i * INPUT_HEIGHT * INPUT_WIDTH + j * INPUT_WIDTH + k] > 0) ? 1 : 0;
            }
        }
    }

    return images;
}

std::vector<int> load_train_labels()
{
    std::ifstream file_train_label("test/mnist_data/train-labels.idx1-ubyte", std::ios::binary);

    if (!file_train_label.is_open()) {
        std::cerr << "Failed to open label file\n";
    }

    // Skip the first 8 bytes of the file, which contain header information
    file_train_label.seekg(8);

    // Read the contents of the file into a vector
    std::vector<unsigned char> train_label((std::istreambuf_iterator<char>(file_train_label)),
                                    std::istreambuf_iterator<char>());

    // Create a 1D array to store the label data
    std::vector<int> labels(NUM_TRAIN_LABELS);

    // Fill the label data from the vector
    for (int i = 0; i < NUM_TRAIN_LABELS; i++) {
        labels[i] = (int)train_label[i];
    }

    return labels;
}


vector<vector<vector<int>>> load_test_data()
{
    ifstream file_test_data("test/mnist_data/t10k-images.idx3-ubyte", std::ios::binary);

    if (!file_test_data.is_open()) {
        std::cerr << "Failed to open image file\n";
    }

    // Skip the first 16 bytes of the file, which contain header information
    file_test_data.seekg(16);

    // Read the contents of the file into a vector
    std::vector<unsigned char> train_data((std::istreambuf_iterator<char>(file_test_data)),
                                    std::istreambuf_iterator<char>());

    // Create a 3D array to store the image data
    std::vector<std::vector<std::vector<int>>> images(NUM_TEST_IMAGES, std::vector<std::vector<int>>(INPUT_HEIGHT, std::vector<int>(INPUT_WIDTH)));

    // Fill the image data from the vector
    for (int i = 0; i < NUM_TEST_IMAGES; i++) {
        for (int j = 0; j < INPUT_HEIGHT; j++) {
            for (int k = 0; k < INPUT_WIDTH; k++) {
                images[i][j][k] = ((int)train_data[i * INPUT_HEIGHT * INPUT_WIDTH + j * INPUT_WIDTH + k] > 0) ? 1 : 0;
            }
        }
    }

    return images;
}

std::vector<int> load_test_labels()
{
    std::ifstream file_test_label("test/mnist_data/t10k-labels.idx1-ubyte", std::ios::binary);

    if (!file_test_label.is_open()) {
        std::cerr << "Failed to open label file\n";
    }

    // Skip the first 8 bytes of the file, which contain header information
    file_test_label.seekg(8);

    // Read the contents of the file into a vector
    std::vector<unsigned char> train_label((std::istreambuf_iterator<char>(file_test_label)),
                                    std::istreambuf_iterator<char>());

    // Create a 1D array to store the label data
    std::vector<int> labels(NUM_TEST_LABELS);

    // Fill the label data from the vector
    for (int i = 0; i < NUM_TEST_LABELS; i++) {
        labels[i] = (int)train_label[i];
    }

    return labels;
}

vector<vector<vector<int>>> load_cnn_weight()
{
    ifstream file("test/model_params/cnn_weight.txt");
    if (!file.is_open()) {
        cerr << "Failed to open the file" << endl;
    }

    int rows, columns, depth;
    file >> rows >> columns >> depth;

    // Create a 3D vector to store the array
    vector<vector<vector<int>>> arr(depth, vector<vector<int>>(rows, vector<int>(columns)));

    // Read the data from the file and store it in the 3D vector
    for (int d = 0; d < depth; ++d) {
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < columns; ++j) {
                file >> arr[d][i][j];
            }
        }
    }

    // Close the file
    file.close();

    return arr;
}

vector<int> load_cnn_bias()
{
    ifstream file("test/model_params/cnn_bias.txt");
    if (!file.is_open()) {
        cerr << "Failed to open the file" << endl;
    }

    int columns;
    file >> columns;

    // Create a 1D vector to store the array
    vector<int> arr(columns);

    // Read the data from the file and store it in the 1D vector
    for (int j = 0; j < columns; ++j) {
        file >> arr[j];
    }

    // Close the file
    file.close();

    return arr;
}

vector<vector<int>> load_fc_weight()
{
    ifstream file("test/model_params/fc_weight.txt");
    if (!file.is_open()) {
        cerr << "Failed to open the file" << endl;
    }

    int rows, columns;
    file >> rows >> columns;

    // Create a 2D vector to store the array
    vector<vector<int>> arr(rows, vector<int>(columns));

    // Read the data from the file and store it in the 2D vector
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
            file >> arr[i][j];
        }
    }

    // Close the file
    file.close();

    return arr;
}

vector<int> load_fc_bias()
{
    ifstream file("test/model_params/fc_bias.txt");
    if (!file.is_open()) {
        cerr << "Failed to open the file" << endl;
    }

    int columns;
    file >> columns;

    // Create a 1D vector to store the array
    vector<int> arr(columns);

    // Read the data from the file and store it in the 1D vector
    for (int j = 0; j < columns; ++j) {
        file >> arr[j];
    }

    // Close the file
    file.close();

    return arr;
}


vector<vector<vector<int>>> conv2d(const vector<vector<int>> &input, const vector<vector<vector<int>>> &kernels, vector<int> &bias) {

    int in_height = input.size();
    int in_width = input[0].size();
    int out_height = (in_height - FILTER_SIZE) / STRIDE + 1;
    int out_width = (in_width - FILTER_SIZE) / STRIDE + 1;
    vector<vector<vector<int>>> outputs(NUM_FILTERS, vector<vector<int>>(out_height, vector<int>(out_width)));

    for (int c = 0; c < NUM_FILTERS; c++) {
        for (int i = 0; i < out_height; i++) {
            for (int j = 0; j < out_width; j++) {
                int sum = 0;
                for (int m = 0; m < FILTER_SIZE; m++) {
                    for (int n = 0; n < FILTER_SIZE; n++) {
                        int x = i * STRIDE + m;
                        int y = j * STRIDE + n;
                        sum += input[x][y] * kernels[c][m][n];
                    }
                }
                outputs[c][i][j] = sum + bias[c];
            }
        }
    }
    return outputs;
}

vector<vector<vector<int>>> max_pooling2d(const vector<vector<vector<int>>> &input) {
  
    int in_height = input[0].size();
    int in_width = input[0][0].size();
    int out_height = in_height / WINDOW_SIZE;
    int out_width = in_width / WINDOW_SIZE;

    vector<vector<vector<int>>> outputs(NUM_FILTERS, vector<vector<int>>(out_height, vector<int>(out_width)));

    for (int c = 0; c < NUM_FILTERS; c++) {
        for (int i = 0; i < out_height; i++) {
            for (int j = 0; j < out_width; j++) {
                int max_val = -1000;
                for (int m = 0; m < WINDOW_SIZE; m++) {
                    for (int n = 0; n < WINDOW_SIZE; n++) {
                        int x = i * WINDOW_SIZE + m;
                        int y = j * WINDOW_SIZE + n;
                        max_val = max(max_val, input[c][x][y]);
                    }
                }
                outputs[c][i][j] = max_val;
            }
        }
    }
    return outputs;
}

vector<int> flatten(const vector<vector<vector<int>>> &inputs) {
  
    vector<int> outputs;

    for (const auto &channel : inputs) {
        for (const auto &row : channel) {
            for (const auto &val : row) {
                outputs.push_back(val);
            }
        }
    }
  return outputs;
}


vector<int> fc_layer(const vector<int>& input, const vector<vector<int>>& weights, const vector<int>& biases) {
    int input_size = input.size();
    int output_size = biases.size();

    vector<int> output(output_size, 0);

    for (int i = 0; i < output_size; i++) {
        for (int j = 0; j < input_size; j++) {
            output[i] += input[j] * weights[j][i];
        }
        output[i] += biases[i];
    }

    return output;
}

vector<vector<int>> transpose_matrix(vector<vector<int>> original_vector){

    vector<vector<int>> transposed_vector((int)original_vector[0].size(), vector<int>((int)original_vector.size()));
    
    for (int i = 0; i < (int)original_vector.size(); ++i) {
        for (int j = 0; j < (int)original_vector[0].size(); ++j) {
            transposed_vector[j][i] = original_vector[i][j];
        }
    }
    
    return transposed_vector;
}

int relu(int x) {
    return x > 0 ? 1 : 0;
}

vector<vector<int>> relu(const vector<vector<int>> &inputs) {
    int width = inputs.size();
    int height = inputs[0].size();

    vector<vector<int>> outputs(width, vector<int>(height, 0));

    for (int i = 0; i < width; i++)
    {
        for (int j=0; j < height; j++)
        {
            outputs[i][j] = relu(inputs[i][j]);
        }
    }

        return outputs;
}

vector<vector<vector<int>>> relu(const vector<vector<vector<int>>> &inputs) {
  int depth = inputs.size();
  int height = inputs[0].size();
  int width = inputs[0][0].size();

  vector<vector<vector<int>>> outputs(depth, vector<vector<int>>(height, vector<int>(width)));

    for (int d = 0; d < depth; d++)
    {
        for (int i = 0; i < height; i++)
        {
            for (int j=0; j < width; j++)
            {
                outputs[d][i][j] = relu(inputs[d][i][j]);
            }
        }
    }
  
    return outputs;
}

int main()
{
    auto test_images = load_test_data();
    auto test_labels = load_test_labels();

    cout << "Data loaded . . ." << endl;

    auto cnn_kernel = load_cnn_weight();
    auto cnn_bias = load_cnn_bias();
    auto fc_weight = load_fc_weight();
    auto fc_weightT = transpose_matrix(fc_weight);
    auto fc_bias = load_fc_bias();

    cout << "Params loaded . . ." << endl;
    
    cout << "Original label: " << test_labels[0] << endl;
    auto output = conv2d(test_images[0], cnn_kernel, cnn_bias);
    cout << "Conv done . . ." << endl;

    output = relu(output);
    cout << "Relu done . . ." << endl;

    output = max_pooling2d(output);
    cout << "Max-pool done . . ." << endl;

    auto fc_in = flatten(output);
    cout << "Flatten done . . ." << endl;

    auto fc_out = fc_layer(fc_in, fc_weightT, fc_bias);
    cout << "FC done . . ." << endl;

    for (int i=0; i < (int)fc_out.size(); i++)
    {
       cout << fc_out[i] << " ";
    }

    cout << "\n\n";

    return 0;
}