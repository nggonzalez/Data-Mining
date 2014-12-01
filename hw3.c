#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>

typedef struct point {
  double x, y;
  int index;
} Point;

typedef struct neighbor {
  int index;
  double distance;
} Neighbor;

typedef struct box {
  Point lowerLeftCorner, upperRightCorner;
  int index,
      parent,
      child1,
      child2,
      child3,
      child4,
      permutationStart,
      permutationEnd,
      searched;
} Box;

/*******************************************************************************
 *******************************************************************************
 *                                METHOD HEADERS                               *
 *******************************************************************************
 ******************************************************************************/
int main(int argc, char *argv[]);
void printBox(Box *b);
void seek(double* a, int n, int k, int* iz);
void buildQtree(Box ***qtree, int *qtreeSize, Point *points, int n, int k,
                int *permutations);
void sort(int *permutations, Point *points, int *indices, int permStart,
          int permEnd, Box child1, Box child2, Box child3, Box child4, int n);
int pointInBox(Point p, Box b);
double generateSearchRadius(Point p, Box *currentBox, Box **qtree);
void findNeighbors(Box **qtree, Box **neighboringLeaves, Point p, double radius,
                   Box *currentBox, int *neighboringLeafIndex, int n);
int overlap(Point p, double radius, Box node);
void findKthDistance(double *distances, int n, int k, double* kth);
void seekHelper(Point* points, int n, int k, int* iz, Point current);
void seek_naive(double* a, int n, int k, int* iz);
void fillPoints(Point* points, double *a, int n);
void findNeighborsNaive(Point* points, Point current, int *sub_iz, int n, int k);
void putNeighborsIntoIz(int* iz, int* sub_iz, int index, int k);
int samePoint(int p1, int p2);
double pointDistance(Point p1, Point p2);
double distance(double x1, double x2, double y1, double y2);
void sortNeighborsByDistance(Neighbor* distances, int n);

int main(int argc, char *argv[]) {
  // int n, k, *iz, *iz_naive, i, j, h;
  // double *a;

  // n = 1000;
  // k = 10;
  // a = malloc(n * 2 * sizeof(double));
  // iz = calloc(n * k, sizeof(int));
  // iz_naive = calloc(n * k, sizeof(int));

  // for(i = 0; i < 2*n; i++){
  //   a[i] = ((double) rand() / (double) RAND_MAX -.5) * 100;
  // }

  // seek_naive(a, n, k, iz_naive);
  // seek(a, n, k, iz);
  // for(i = 0; i < n; i++) {
  //   printf("(%d)[%.5f, %.5f]:\t", i, a[2 * i], a[2 * i + 1]);
  //   for(j = 0; j < k; j++){
  //     printf("%d\t", iz[i * k + j]);
  //   }
  //   printf("||\t");
  //   for(j = 0; j < k; j++){
  //     printf("%d\t", iz_naive[i * k + j]);
  //   }
  //   printf("\n");
  // }
  return 1;
}

/*******************************************************************************
 *******************************************************************************
 *                                SEEK FAST CODE                               *
 *******************************************************************************
 ******************************************************************************/
void seek(double* a, int n, int k, int* iz) {
  int i, j, h;

  Point* points = malloc(n * sizeof(Point));
  fillPoints(points, a, n);

  int *permutations = malloc(n * sizeof(int));
  for(i = 0; i < n; i++) {
    permutations[i] = i;
  }

  int *qtreeSize = malloc(sizeof(int));
  *qtreeSize = n;
  double radius;
  Box **qtree, *currentBox;
  qtree = calloc((*qtreeSize), sizeof(Box*));

  buildQtree(&qtree, qtreeSize, points, n, k, permutations);

  Box **searchQueue = malloc((*qtreeSize) * sizeof(Box *));
  Box **leaves = malloc((*qtreeSize) * sizeof(Box *));
  int searchIndex = 1, leafCount = 0;
  searchQueue[0] = qtree[0];

  for(i = 0; i < searchIndex; i++) {
    Box *node = searchQueue[i];
    if(node->child1 > -1) {
      searchQueue[searchIndex] = qtree[node->child1];
      searchQueue[searchIndex + 1] = qtree[node->child2];
      searchQueue[searchIndex + 2] = qtree[node->child3];
      searchQueue[searchIndex + 3] = qtree[node->child4];
      searchIndex += 4;
    } else { //leaf
      leaves[leafCount++] = node;
    }
  }

  // Store what leaf box the point is in at the end of building the tree
  Box **pointsToLeaves = malloc(n * sizeof(Box *));
  for(i = 0; i < leafCount; i++) {
    for(j = leaves[i]->permutationStart; j < leaves[i]->permutationEnd; j++) {
      pointsToLeaves[permutations[j]] = leaves[i];
    }
  }

  for(i = 0; i < n; i++) {
    currentBox = pointsToLeaves[i];
    double radius = generateSearchRadius(points[i], currentBox, qtree);
    Box **neighboringLeaves = calloc(leafCount, sizeof (Box *));
    int neighboringLeafCount = 0;
    int *neighboringLeafIndex = malloc(sizeof(int));
    *neighboringLeafIndex = 0;

    findNeighbors(qtree, neighboringLeaves, points[i], radius, currentBox,
                  neighboringLeafIndex, n);

    int subN = 0;
    Point *subPoints = malloc(n * sizeof(Point));
    for(j = 0; j < (*neighboringLeafIndex); j++) {
      for(h = neighboringLeaves[j]->permutationStart; h < neighboringLeaves[j]->permutationEnd; h++) {
        subPoints[subN++] = points[permutations[h]];
      }
    }

    int *sub_iz = malloc(k * sizeof(int));
    // findNeighborsNaive(subPoints, points[i], sub_iz, subN, k);
    seekHelper(subPoints, subN, k, sub_iz, points[i]);
    putNeighborsIntoIz(iz, sub_iz, i, k);
  }
}

void buildQtree(Box ***qtree, int *qtreeSize, Point *points, int n, int k,
                int *permutations) {
  double xMax = INT_MIN, yMax = INT_MIN;
  double xMin = INT_MAX, yMin = INT_MAX;

  int i;
  for(i = 0; i < n; i++) {
    Point p = points[i];
    if(p.x > xMax) {
      xMax = p.x;
    }
    if(p.x < xMin) {
      xMin = p.x;
    }
    if(p.y > yMax) {
      yMax = p.y;
    }
    if(p.y < yMin) {
      yMin = p.y;
    }
  }

  // root square
  Box *root = malloc(sizeof(Box));
  root->upperRightCorner.x = xMax;
  root->upperRightCorner.y = yMax;
  root->lowerLeftCorner.x = xMin;
  root->lowerLeftCorner.y = yMin;
  // Initial permutation matrix - points in the box
  root->permutationStart = 0;
  root->permutationEnd = n;

  // Parents
  root->parent = 0;
  root->searched = 0;
  root->index = 0;

  // Set root to first box in array for the quad tree
  (*qtree)[0] = root;

  int currentBoxIndex = 0, maxBoxIndex = 0, minBoxIndex = 0;
  Box *currentBox;

  while(currentBoxIndex >= 0) {
    // Make tree bigger
    if(maxBoxIndex > (*qtreeSize) - 10) {
      (*qtreeSize) *= 2;
      Box **temp = calloc((*qtreeSize), sizeof(Box*));
      for(i = 0; i < maxBoxIndex + 4; i++) {
        temp[i] = (*qtree)[i];
      }
      free(*qtree);
      *qtree = temp;
    }

    currentBoxIndex = -1;

    // Unsearched boxes?
    for(i = minBoxIndex; i < maxBoxIndex + 1; i++) {
      minBoxIndex = i;
      if((*qtree)[i]->searched == 0) {
        currentBoxIndex = i;
        currentBox = (*qtree)[i];
        break;
      }
    }

    // If more than k points in the box, then give the box chilun
    if(currentBox->permutationEnd - currentBox->permutationStart > k) {
      Box *child1, *child2, *child3, *child4;
      child1 = malloc(sizeof(Box));
      child2 = malloc(sizeof(Box));
      child3 = malloc(sizeof(Box));
      child4 = malloc(sizeof(Box));

      // Set new chiluns chilun to "undefined"
      child1->child1 = child1->child2 = child1->child3 = child1->child4 = -1;
      child2->child1 = child2->child2 = child2->child3 = child2->child4 = -1;
      child3->child1 = child3->child2 = child3->child3 = child3->child4 = -1;
      child4->child1 = child4->child2 = child4->child3 = child4->child4 = -1;

      // Initialize search value
      child1->searched = child2->searched = child3->searched = child4->searched = 0;

      // Initialize parent value
      child1->parent = child2->parent = child3->parent = child4->parent = currentBoxIndex;

      // Set chiluns' boundaries
      child1->lowerLeftCorner.x = currentBox->lowerLeftCorner.x;
      child1->lowerLeftCorner.y = currentBox->lowerLeftCorner.y+((currentBox->upperRightCorner.y - currentBox->lowerLeftCorner.y)/2);
      child1->upperRightCorner.x = currentBox->lowerLeftCorner.x+((currentBox->upperRightCorner.x - currentBox->lowerLeftCorner.x)/2);
      child1->upperRightCorner.y = currentBox->upperRightCorner.y;

      child2->lowerLeftCorner.x = currentBox->lowerLeftCorner.x+((currentBox->upperRightCorner.x - currentBox->lowerLeftCorner.x)/2);
      child2->lowerLeftCorner.y = currentBox->lowerLeftCorner.y+((currentBox->upperRightCorner.y - currentBox->lowerLeftCorner.y)/2);
      child2->upperRightCorner.x = currentBox->upperRightCorner.x;
      child2->upperRightCorner.y = currentBox->upperRightCorner.y;

      child3->lowerLeftCorner.x = currentBox->lowerLeftCorner.x;
      child3->lowerLeftCorner.y = currentBox->lowerLeftCorner.y;
      child3->upperRightCorner.x = currentBox->lowerLeftCorner.x+((currentBox->upperRightCorner.x - currentBox->lowerLeftCorner.x)/2);
      child3->upperRightCorner.y = currentBox->lowerLeftCorner.y+((currentBox->upperRightCorner.y - currentBox->lowerLeftCorner.y)/2);

      child4->lowerLeftCorner.x = currentBox->lowerLeftCorner.x+((currentBox->upperRightCorner.x - currentBox->lowerLeftCorner.x)/2);
      child4->lowerLeftCorner.y = currentBox->lowerLeftCorner.y;
      child4->upperRightCorner.x = currentBox->upperRightCorner.x;
      child4->upperRightCorner.y = currentBox->lowerLeftCorner.y+((currentBox->upperRightCorner.y - currentBox->lowerLeftCorner.y)/2);

      int *indices = malloc(5 * sizeof(int));

      sort(permutations, points, indices, currentBox->permutationStart, currentBox->permutationEnd, *child1, *child2, *child3, *child4, n);

      // Assign indices to children
      child1->permutationStart = indices[0];
      child1->permutationEnd = indices[1];
      child2->permutationStart = indices[1];
      child2->permutationEnd = indices[2];
      child3->permutationStart = indices[2];
      child3->permutationEnd = indices[3];
      child4->permutationStart = indices[3];
      child4->permutationEnd = indices[4];

      // Add children to the qtree
      (*qtree)[maxBoxIndex + 1] = child1;
      (*qtree)[maxBoxIndex + 2] = child2;
      (*qtree)[maxBoxIndex + 3] = child3;
      (*qtree)[maxBoxIndex + 4] = child4;

      currentBox->child1 = child1->index = maxBoxIndex + 1;
      currentBox->child2 = child2->index = maxBoxIndex + 2;
      currentBox->child3 = child3->index = maxBoxIndex + 3;
      currentBox->child4 = child4->index = maxBoxIndex + 4;

      maxBoxIndex += 4;
    }

    // Checked box and gave it children if necessary
    if(currentBoxIndex >= 0) {
      currentBox->searched = 1;
    }
  }
}

void sort(int *permutations, Point *points, int *indices, int permStart,
          int permEnd, Box child1, Box child2, Box child3, Box child4, int n) {
  int child1Index, child2Index, child3Index, child4Index;
  child1Index = child2Index = child3Index = child4Index = 0;

  int *child1Array = malloc(n * sizeof(int));
  int *child2Array = malloc(n * sizeof(int));
  int *child3Array = malloc(n * sizeof(int));
  int *child4Array = malloc(n * sizeof(int));

  int i;
  for(i = permStart; i < permEnd; i++) {
    if(pointInBox(points[permutations[i]], child1)) {
      child1Array[child1Index++] = permutations[i];
    } else if(pointInBox(points[permutations[i]], child2)) {
      child2Array[child2Index++] = permutations[i];
    } else if(pointInBox(points[permutations[i]], child3)) {
      child3Array[child3Index++] = permutations[i];
    } else if(pointInBox(points[permutations[i]], child4)) {
      child4Array[child4Index++] = permutations[i];
    }
  }

  indices[0] = i = permStart;
  int indexSum = permStart;
  if(child1Index > 0) {
    indexSum += child1Index;
    for(; i < indexSum; i++) {
      permutations[i] = child1Array[i - permStart];
    }
  }

  indices[1] = i;
  if(child2Index > 0) {
    indexSum += child2Index;
    for(; i < indexSum; i++) {
      permutations[i] = child2Array[i - child1Index - permStart];
    }
  }

  indices[2] = i;
  if(child3Index > 0) {
    indexSum += child3Index;
    for(; i < indexSum; i++) {
      permutations[i] = child3Array[i - child1Index - child2Index - permStart];
    }
  }

  indices[3] = i;
  if(child4Index > 0) {
    indexSum += child4Index;
    for(; i < indexSum; i++) {
      permutations[i] = child4Array[i - child1Index - child2Index - child3Index - permStart];
    }
  }

  indices[4] = indexSum;
}

int pointInBox(Point p, Box b) {
  return p.x >= b.lowerLeftCorner.x && p.y >= b.lowerLeftCorner.y &&
         p.x <= b.upperRightCorner.x && p.y <= b.upperRightCorner.y;
}

double generateSearchRadius(Point p, Box *currentBox, Box **qtree) {
  Point upperLeftCorner, lowerRightCorner;
  Box *parent = qtree[currentBox->parent];

  upperLeftCorner.x = parent->lowerLeftCorner.x;
  upperLeftCorner.y = parent->upperRightCorner.y;
  lowerRightCorner.x = parent->upperRightCorner.x;
  lowerRightCorner.y = parent->lowerLeftCorner.y;

  return fmax(pointDistance(p, upperLeftCorner),
              fmax(pointDistance(p, parent->upperRightCorner),
                   fmax(pointDistance(p, lowerRightCorner),
                        pointDistance(p, parent->lowerLeftCorner))));
}

void findNeighbors(Box **qtree, Box **neighboringLeaves, Point p, double radius,
                   Box *currentBox, int *neighboringLeafIndex, int n) {
  Box **searchQueue = malloc(n * sizeof(Box *));
  int searchIndex = 1;
  searchQueue[0] = qtree[0];

  int i;
  for(i = 0; i < searchIndex; i++) {
    Box *node = searchQueue[i];
    if(!overlap(p, radius, *(node))) {
      continue;
    }
    if(node->child1 != -1) {
      searchQueue[searchIndex++] = qtree[node->child1];
      searchQueue[searchIndex++] = qtree[node->child2];
      searchQueue[searchIndex++] = qtree[node->child3];
      searchQueue[searchIndex++] = qtree[node->child4];
    } else { //leaf
      neighboringLeaves[(*neighboringLeafIndex)++] = node;
    }
  }
}

int overlap(Point p, double radius, Box node) {
  if(pointInBox(p, node)) {
    return 1;
  } else {
    double xMidCoord = (node.upperRightCorner.x + node.lowerLeftCorner.x) / 2;
    double yMidCoord = (node.upperRightCorner.y + node.lowerLeftCorner.y) / 2;
    double width = node.upperRightCorner.x - node.lowerLeftCorner.x;
    double height = node.upperRightCorner.y - node.lowerLeftCorner.y;

    Point circle;
    circle.x = fabs(p.x - xMidCoord);
    circle.y = fabs(p.y - yMidCoord);

    if (circle.x > (width / 2 + radius)) {
      return 0;
    } else if (circle.y > (height / 2 + radius)) {
      return 0;
    } else if (circle.x <= (width / 2)) {
      return 1;
    } else if (circle.y <= (height / 2)) {
      return 1;
    }

    double cornerSquare = pow((circle.x - width/2), 2) + pow((circle.y - height/2), 2);
    return (cornerSquare <= pow(radius, 2));
  }
}


/*******************************************************************************
 *******************************************************************************
 *                               SEEK NAIVE CODE                               *
 *******************************************************************************
 ******************************************************************************/

void seek_naive(double* a, int n, int k, int* iz) {
  Point* points = malloc(n * sizeof(Point));
  fillPoints(points, a, n);

  int i;

  for(i = 0; i < n; i++) {
    int *sub_iz = malloc(k * sizeof(int));
    // findNeighborsNaive(points, points[i], sub_iz, n, k);
    seekHelper(points, n, k, sub_iz, points[i]);
    putNeighborsIntoIz(iz, sub_iz, i, k);
    free(sub_iz);
  }
}

void findNeighborsNaive(Point* points, Point current, int *sub_iz, int n, int k) {
  Neighbor* distances = malloc((n - 1) * sizeof(Neighbor));

  int i, index;
  for (index = i = 0; index < n - 1 && i < n; index++) {
    if(samePoint(points[i].index, current.index)) {
      i++;
    }
    distances[index].index = points[i].index;
    distances[index].distance = pointDistance(points[i], current);
    i++;
  }

  sortNeighborsByDistance(distances, n -  1);
  for(i = 0; i < k; i++) {
    sub_iz[i] = distances[i].index;
  }

  free(distances);
}

void sortNeighborsByDistance(Neighbor* distances, int n) {
  int k, i;

  for(k = n - 1; k > 0; k--) {
    for(i = 0; i < k; i++) {
      if(distances[i].distance > distances[i + 1].distance) {
        Neighbor swap = distances[i];
        distances[i] = distances[i + 1];
        distances[i + 1] = swap;
      }
    }
  }
}

double pointDistance(Point p1, Point p2) {
  return distance(p1.x, p2.x, p1.y, p2.y);
}

double distance(double x1, double x2, double y1, double y2) {
  return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
}

int samePoint(int p1, int p2) {
  return p1 == p2;
}

void fillPoints(Point* points, double *a, int n) {
  int i, index = 0;
  for(i = 0; i < 2 * n; i += 2) {
    points[index].x = a[i];
    points[index].y = a[i + 1];
    points[index++].index = i / 2;
    // printf("%f %f, index: %d\n", points[index - 1].x, points[index - 1].y, points[index - 1].index);
  }
}

void putNeighborsIntoIz(int* iz, int* sub_iz, int index, int k) {
  int i;
  for(i = 0; i < k; i++) {
    iz[k * index + i] = sub_iz[i]; // Move past other points neighbors
  }
}

void printBox(Box *b){
  if (b->permutationStart >= b->permutationEnd) {
    return;
  }
  printf("x values: [%f, %f]\n", b->lowerLeftCorner.x, b->upperRightCorner.x);
  printf("y values: [%f, %f]\n", b->lowerLeftCorner.y, b->upperRightCorner.y);
  printf("points: %d - %d\n", b->permutationStart+1, b->permutationEnd);
  printf("check: %d\n", b->searched);
  printf("parent: %d\n", b->parent);
  printf("index: %d\n", b->index);
  printf("children: [%d, %d, %d, %d]\n", b->child1, b->child2, b->child3, b->child4);
  printf("\n");
}

void findKthDistance(double *distances, int n, int k, double* kth) {
  int i, start;
  double temp;

  for (start = i = 0; i < n - 1; i++) {
    if (distances[i] > distances[n - 1]) {
      continue;
    }
    temp = distances[i];
    distances[i] = distances[start];
    distances[start] = temp;
    start++;
  }
  temp = distances[n - 1];
  distances[n - 1] = distances[start];
  distances[start] = temp;

  if (k == start) {
      *(kth) = distances[start];
      return;
  } else if(start > k) {
    findKthDistance(distances, start, k, kth);
  } else {
    findKthDistance(distances + start, n - start, k - start, kth);
  }
}

void seekHelper(Point* points, int n, int k, int* iz, Point current) {
  Point* neighbors = malloc(k * sizeof(Point));
  double* distances = malloc(n * sizeof(double));
  double* selection = malloc(n * sizeof(double));

  int i;
  for (i = 0; i < n; i++) {
    distances[i] = pointDistance(points[i], current);
    selection[i] = distances[i];
  }

  double* kth = malloc(sizeof(double));
  findKthDistance(selection, n, k, kth);

  int neighborIndex = 0;
  for (i = 0; i < n; i++) {
    if (samePoint(points[i].index, current.index)) {
      continue;
    } else if (distances[i] <= *(kth)) {
      neighbors[neighborIndex++] = points[i];
    }
  }

  for(i = 0; i < k; i++) {
    iz[i] = neighbors[i].index;
  }

  return;
}