import edu.princeton.cs.algs4.Picture;
import java.awt.Color;
import java.util.Arrays;

public class SeamCarver {
    private Picture myPicture;
    private Double[][] energyMatrix;
    private Double[][] rowDistTo;
    private Double[][] colDistTo;
    private int[][] rowEdgeTo;
    private int[][] columnEdgeTo;
    private Color[][] pictureMatrix;
    private boolean isTransposed;

    // create a seam carver object based on the given picture
    public SeamCarver(Picture picture) {
        // first validate the argument passed to the constructor
        validateConstructor(picture);
        myPicture = new Picture(picture);
        int pictureHeight = myPicture.height();
        int pictureWidth = myPicture.width();
        // java matrix representation is r x c
        pictureMatrix = new Color[pictureHeight][pictureWidth];

        // build a picture matrix of color objects
        for (int i = 0; i < pictureHeight; i++) {
            for (int j = 0; j < pictureWidth; j++) {
                pictureMatrix[i][j] = myPicture.get(j, i);
                // int r = (myPicture.getRGB(i, j) >> 16) & 0xFF;
                // int g = (myPicture.getRGB(i, j) >> 8) & 0xFF;
                // int b = (myPicture.getRGB(i, j) >> 0) & 0xFF;
            }
        }
        // r x c ==> h x w
        energyMatrix = new Double[pictureHeight][pictureWidth];

        calculateTotalPixelEnergy(pictureWidth, pictureHeight);

        // initialize distTo arrays
        initializeDistToArray(pictureWidth, pictureHeight);

        // initialize edgeTo arrays
        initializeEdgeToArrays(pictureWidth, pictureHeight);
    }

    private void calculateTotalPixelEnergy(int pWidth, int pHeight) {
        // calculate the energy of each pixel
        for (int c = 0; c < pWidth; c++) {
            for (int r = 0; r < pHeight; r++) {
                calculatePixelEnergy(c, r);
            }
        }
    }

    private void calculatePixelEnergy(int c, int r) {
        // the energy of a border pixel is 1000
        // java matrix is r x c
        // but pixel representation is c x r
        if (c == 0) {
            energyMatrix[r][c] = 1000.0;
        }
        else if (c >= width() - 1) {
            energyMatrix[r][c] = 1000.0;
        }
        else if (r == 0) {
            energyMatrix[r][c] = 1000.0;
        }
        else if (r >= height() - 1) {
            energyMatrix[r][c] = 1000.0;
        }
        else {
            Color leftPixel = pictureMatrix[r][c - 1];
            Color rightPixel = pictureMatrix[r][c + 1];
            Color topPixel = pictureMatrix[r - 1][c];
            Color bottomPixel = pictureMatrix[r + 1][c];
            // calculate x central difference in all RGB colors
            int xRed = leftPixel.getRed() - rightPixel.getRed();
            int xGreen = leftPixel.getGreen() - rightPixel.getGreen();
            int xBlue = leftPixel.getBlue() - rightPixel.getBlue();
            // calculate x-gradient square
            int xGradientSquare = (xRed * xRed) + (xGreen * xGreen) + (xBlue * xBlue);

            // calculate y central difference in all RGB colors
            int yRed = topPixel.getRed() - bottomPixel.getRed();
            int yGreen = topPixel.getGreen() - bottomPixel.getGreen();
            int yBlue = topPixel.getBlue() - bottomPixel.getBlue();
            // calculate y-gradient square
            int yGradientSquare = (yRed * yRed) + (yGreen * yGreen) + (yBlue * yBlue);

            energyMatrix[r][c] = Math.sqrt(xGradientSquare + yGradientSquare);
        }
    }

    private void initializeDistToArray(int pWidth, int pHeight) {
        rowDistTo = new Double[pHeight][pWidth];
        // initialize the distTo array
        for (int i = 0; i < pWidth; i++) {
            for (int j = 0; j < pHeight; j++) {
                // the distance to the top row is 0
                if (j == 0) {
                    rowDistTo[j][i] = 0.0;
                }
                else {
                    rowDistTo[j][i] = Double.POSITIVE_INFINITY;
                }
            }
        }

        colDistTo = new Double[pHeight][pWidth];
        // initialize the distTo array
        for (int i = 0; i < pWidth; i++) {
            for (int j = 0; j < pHeight; j++) {
                // the energy of a border pixel is 1000
                if (i == 0) {
                    colDistTo[j][i] = 0.0;
                }
                else {
                    colDistTo[j][i] = Double.POSITIVE_INFINITY;
                }
            }
        }
    }

    private void initializeEdgeToArrays(int pWidth, int pHeight) {
        rowEdgeTo = new int[pHeight][pWidth];
        columnEdgeTo = new int[pHeight][pWidth];
    }

    // current picture
    public Picture picture() {
        // if the current picture matrix is a transposed matrix
        // then transpose it back to the originel
        if (isTransposed) {
            pictureMatrix = castToColor(transpose(pictureMatrix));
            energyMatrix = castToDouble(transpose(energyMatrix));
        }

        // update the picture from the pictureMatrix
        myPicture = new Picture(pictureMatrix[0].length, pictureMatrix.length);
        for (int i = 0; i < pictureMatrix[0].length; i++) {
            for (int j = 0; j < pictureMatrix.length; j++) {
                myPicture.set(i, j, pictureMatrix[j][i]);
            }
        }
        return myPicture;
    }

    // width of current picture
    public int width() {
        return pictureMatrix[0].length;
    }

    // height of current picture
    public int height() {
        return pictureMatrix.length;
    }

    // energy of pixel at column x and row y
    public double energy(int x, int y) {
        // COMPLETED: I need to first check whether it is a transposed matrix or not.
        // if the current picture matrix is a transposed matrix
        // then transpose it back to the originel
        if (isTransposed) {
            pictureMatrix = castToColor(transpose(pictureMatrix));
            energyMatrix = castToDouble(transpose(energyMatrix));
        }
        if (x < 0 || x > width() - 1 || y < 0 || y > height() - 1) {
            throw new IllegalArgumentException("An index is outside of its prescribed range");
        }
        return energyMatrix[y][x];
    }

    // sequence of indices for horizontal seam
    public int[] findHorizontalSeam() {
        // if the pictureMatrix has been transposed earlier, then transpose it back to the original matrix
        if (isTransposed) {
            pictureMatrix = castToColor(transpose(pictureMatrix));
            energyMatrix = castToDouble(transpose(energyMatrix));
            isTransposed = !isTransposed;
        }

        for (int c = 0; c < width(); c++) {
            for (int r = 0; r < height(); r++) {
                // relax all the edges from vertices in column 0
                // for vertex 0,1, the edges are 1,0 1,1 1,2
                relaxColumnEdge(c, r, energyMatrix[r][c]);
            }

        }

        // go through all of the vertices(rows) in the last column and pick the one with the shortest distance
        // Completed: how to select one of three shortest distance vertices in the last column
        double shortestDist = Double.POSITIVE_INFINITY;
        int row = -1;
        for (int r = 0; r < height(); r++) {
            if (colDistTo[r][width() - 1] < shortestDist) {
                shortestDist = colDistTo[r][width() - 1];
                row = r;
            }
        }
        // track back recording the rows in the seam
        int[] solution = new int[width()];

        // starting from the last column and working the way back
        int columnIndex = width() - 1;
        while (colDistTo[row][columnIndex] != 0.0) {
            solution[columnIndex] = row;
            row = columnEdgeTo[row][columnIndex];
            columnIndex = columnIndex - 1;
        }
        solution[columnIndex] = row;

        return solution;
    }

    private void relaxColumnEdge(int c, int r, double weight) {
        // if we are at the last column, just return
        if (c == width() - 1) {
            return;
        }
        // if there is only one row (or the height of the image == 1)
        if (r == 0 && r == height() - 1) {
            relaxColumnEdgeRightMiddle(c, r, weight);
        }
        // if we are at the top row, only relax for rightmiddle and rightbottom
        else if (r == 0) {
            relaxColumnEdgeRightMiddle(c, r, weight);
            relaxColumnEdgeRightBottom(c, r, weight);
        }
        // if we are at the bottom row, only relax righttop and rightmiddle
        else if (r == height() - 1) {
            relaxColumnEdgeRightTop(c, r, weight);
            relaxColumnEdgeRightMiddle(c, r, weight);
        }
        // otherwise relax all the vertices to the right (righttop, rightmiddle and rightbottom)
        else {
            relaxColumnEdgeRightTop(c, r, weight);
            relaxColumnEdgeRightMiddle(c, r, weight);
            relaxColumnEdgeRightBottom(c, r, weight);
        }
    }

    private void relaxColumnEdgeRightTop(int c, int r, double weight) {
        if (colDistTo[r - 1][c + 1] > colDistTo[r][c] + weight) {
            colDistTo[r - 1][c + 1] = colDistTo[r][c] + weight;
            columnEdgeTo[r - 1][c + 1] = r;
        }
    }

    private void relaxColumnEdgeRightMiddle(int c, int r, double weight) {
        if (colDistTo[r][c + 1] > colDistTo[r][c] + weight) {
            colDistTo[r][c + 1] = colDistTo[r][c] + weight;
            columnEdgeTo[r][c + 1] = r;
        }
    }

    private void relaxColumnEdgeRightBottom(int c, int r, double weight) {
        if (colDistTo[r + 1][c + 1] > colDistTo[r][c] + weight) {
            colDistTo[r + 1][c + 1] = colDistTo[r][c] + weight;
            columnEdgeTo[r + 1][c + 1] = r;
        }
    }

    // sequence of indices for vertical seam
    public int[] findVerticalSeam() {
        // if the pictureMatrix has been transposed earlier, then transpose it back to the original matrix
        if (isTransposed) {
            pictureMatrix = castToColor(transpose(pictureMatrix));
            energyMatrix = castToDouble(transpose(energyMatrix));
            isTransposed = !isTransposed;
        }

        for (int r = 0; r < height(); r++) {
            for (int c = 0; c < width(); c++) {
                relaxRowEdge(c, r, energyMatrix[r][c]);
            }
        }
        // go through all of the vertices(columns) in the last row and pick the one with the shortest distance
        // TODO: how to select one of three shortest distance vertices in the last row
        double shortestDist = Double.POSITIVE_INFINITY;
        int col = -1;
        for (int c = 0; c < width(); c++) {
            if (rowDistTo[height() - 1][c] < shortestDist) {
                shortestDist = rowDistTo[height() - 1][c];
                col = c;
            }
        }
        // track back recording the columns in the seam
        int[] solution = new int[height()];

        // starting from the last row and working the way back up
        int rowIndex = height() - 1;
        while (rowDistTo[rowIndex][col] != 0.0) {
            solution[rowIndex] = col;
            col = rowEdgeTo[rowIndex][col];
            rowIndex = rowIndex - 1;
        }
        solution[rowIndex] = col;

        return solution;
    }

    private void relaxRowEdge(int c, int r, double weight) {
        // if we are at the last row, just return
        if (r == height() - 1) {
            return;
        }
        // if there is only one column
        if (c == 0 && c == width() - 1) {
            relaxRowEdgeBottomMiddle(c, r, weight);
        }
        // if we are at the leftmost column, only relax for bottommiddle and bottomright
        else if (c == 0) {
            relaxRowEdgeBottomMiddle(c, r, weight);
            relaxRowEdgeBottomRight(c, r, weight);
        }
        // if we are at the rightmost column, only relax bottommiddle and bottomleft
        else if (c == width() - 1) {
            relaxRowEdgeBottomLeft(c, r, weight);
            relaxRowEdgeBottomMiddle(c, r, weight);
        }
        // otherwise relax all the vertices to the bottom (bottomleft, bottommiddle and bottomright)
        else {
            relaxRowEdgeBottomLeft(c, r, weight);
            relaxRowEdgeBottomMiddle(c, r, weight);
            relaxRowEdgeBottomRight(c, r, weight);
        }
    }

    private void relaxRowEdgeBottomRight(int c, int r, double weight) {
        if (rowDistTo[r + 1][c + 1] > rowDistTo[r][c] + weight) {
            rowDistTo[r + 1][c + 1] = rowDistTo[r][c] + weight;
            rowEdgeTo[r + 1][c + 1] = c;
        }
    }

    private void relaxRowEdgeBottomMiddle(int c, int r, double weight) {
        if (rowDistTo[r + 1][c] > rowDistTo[r][c] + weight) {
            rowDistTo[r + 1][c] = rowDistTo[r][c] + weight;
            rowEdgeTo[r + 1][c] = c;
        }
    }

    private void relaxRowEdgeBottomLeft(int c, int r, double weight) {
        if (rowDistTo[r + 1][c - 1] > rowDistTo[r][c] + weight) {
            rowDistTo[r + 1][c - 1] = rowDistTo[r][c] + weight;
            rowEdgeTo[r + 1][c - 1] = c;
        }
    }

    // remove horizontal seam from current picture
    public void removeHorizontalSeam(int[] seam) {

        // if the picture matrix is not transposed, then transpose it
        // transpose the energy matrix too
        if (!isTransposed) {
            pictureMatrix = castToColor(transpose(pictureMatrix));
            energyMatrix = castToDouble(transpose(energyMatrix));
            isTransposed = !isTransposed;
        }
        validateVertex(seam);

        if (pictureMatrix[0].length <= 1) {
            throw new IllegalArgumentException(
                    "The height of the picture is too small to remove a horizontal seam");
        }

        // create an appropriate size dest picture matrix
        Color[][] destPictureMatrix = new Color[pictureMatrix.length][pictureMatrix[0].length - 1];
        Double[][] destEnergyMatrix = new Double[pictureMatrix.length][pictureMatrix[0].length - 1];

        // remove the appropriate indices
        pictureMatrix = (Color[][]) removeSeam(pictureMatrix, destPictureMatrix, seam);
        energyMatrix = (Double[][]) removeSeam(energyMatrix, destEnergyMatrix, seam);

        // transpose back the pictureMatrix, when the picture is needed,
        // otherwise, we might need to do further horizontal seam removal.
        // In which case, there would be no need to transpose back yet
        // update the picture
        // Completed: update the energy matrix, for the pixels that now occupy where seam was removed.
        for (int r = 0; r < pictureMatrix.length; r++) {
            if (seam[r] == pictureMatrix[0].length) {
                calculatePixelEnergy(seam[r] - 1, r);
            }
            else {
                calculatePixelEnergy(seam[r], r);
            }
        }
        initializeDistToArray(pictureMatrix.length, pictureMatrix[0].length);
        initializeEdgeToArrays(pictureMatrix.length, pictureMatrix[0].length);
    }

    private Object[][] transpose(Object[][] matrix) {
        Object[][] transposedPictureMatrix = new Object[matrix[0].length][matrix.length];
        for (int i = 0; i < matrix[0].length; i++) {
            for (int j = 0; j < matrix.length; j++) {
                transposedPictureMatrix[i][j] = matrix[j][i];
            }
        }
        return transposedPictureMatrix;
    }

    private Color[][] castToColor(Object[][] objectMatrixTypeColor) {
        Color[][] matrix = new Color[objectMatrixTypeColor.length][objectMatrixTypeColor[0].length];
        for (int i = 0; i < objectMatrixTypeColor.length; i++) {
            for (int j = 0; j < objectMatrixTypeColor[0].length; j++) {
                matrix[i][j] = (Color) objectMatrixTypeColor[i][j];
            }
        }
        return matrix;
    }

    private Double[][] castToDouble(Object[][] objectMatrixTypeDouble) {
        Double[][] matrix
                = new Double[objectMatrixTypeDouble.length][objectMatrixTypeDouble[0].length];
        for (int i = 0; i < objectMatrixTypeDouble.length; i++) {
            for (int j = 0; j < objectMatrixTypeDouble[0].length; j++) {
                matrix[i][j] = (Double) objectMatrixTypeDouble[i][j];
            }
        }
        return matrix;
    }

    // remove vertical seam from current picture
    public void removeVerticalSeam(int[] seam) {
        if (isTransposed) {
            pictureMatrix = castToColor(transpose(pictureMatrix));
            energyMatrix = castToDouble(transpose(energyMatrix));
            isTransposed = !isTransposed;
        }
        validateVertex(seam);

        if (pictureMatrix[0].length <= 1) {
            throw new IllegalArgumentException(
                    "The width of the picture is too small to remove a vertical seam");
        }
        // new color[width][height]
        Color[][] destPictureMatrix = new Color[pictureMatrix.length][pictureMatrix[0].length - 1];
        Double[][] destEnergyMatrix = new Double[pictureMatrix.length][pictureMatrix[0].length - 1];

        pictureMatrix = (Color[][]) removeSeam(pictureMatrix, destPictureMatrix, seam);
        energyMatrix = (Double[][]) removeSeam(energyMatrix, destEnergyMatrix, seam);

        // Completed: update the energy matrix, for the pixels that now occupy where seam was removed.
        for (int r = 0; r < pictureMatrix.length; r++) {
            if (seam[r] == pictureMatrix[0].length) {
                calculatePixelEnergy(seam[r] - 1, r);
            }
            else {
                calculatePixelEnergy(seam[r], r);
            }
        }

        initializeDistToArray(pictureMatrix[0].length, pictureMatrix.length);
        initializeEdgeToArrays(pictureMatrix[0].length, pictureMatrix.length);
    }

    private Object[][] removeSeam(Object[][] origMatrix, Object[][] destMatrix, int[] seam) {
        for (int r = 0; r < origMatrix.length; r++) {
            int indexToRemove = seam[r];
            System.arraycopy(origMatrix[r], 0, destMatrix[r], 0, indexToRemove);
            System.arraycopy(origMatrix[r], indexToRemove + 1, destMatrix[r],
                             indexToRemove,
                             origMatrix[0].length - (indexToRemove + 1));
        }
        return destMatrix;
    }

    private void validateConstructor(Picture picture) {
        if (picture == null) {
            throw new IllegalArgumentException("An argument to the constructor is null");
        }
    }

    private void validateVertex(int[] seam) {
        if (seam == null) {
            throw new IllegalArgumentException("An argument to the method is null");
        }

        if (seam.length != height()) {
            throw new IllegalArgumentException(
                    "An array argument to the method is of wrong length and does "
                            + "not match the corresponding direction for the picture matrix");
        }

        // throw an exception if an entry is outside of its prescribed range
        for (int i = 0; i < seam.length; i++) {
            // for every row in the pictureMatrix, if the column to be removed is
            // outside of the range of the pictureMatrix column
            if (seam[i] < 0 || seam[i] > width() - 1) {
                throw new IllegalArgumentException(
                        "Invalid seam - an entry to the seam is outside of its prescribed range");
            }
        }


        // throw null if adjacent entries differ by more than one
        for (int i = 0; i < seam.length - 1; i++) {
            if (Math.abs(seam[i] - seam[i + 1]) > 1) {
                throw new IllegalArgumentException(
                        "Invalid seam - adjacent entries differ by more than one");
            }
        }
    }


    //  unit testing (optional)
    public static void main(String[] args) {
        Picture picture = new Picture(args[0]);
        SeamCarver sc = new SeamCarver(picture);

        Double[][] matrix = new Double[3][4];
        Double[][] destMatrix = new Double[3][3];
        matrix[0][0] = 0.0;
        matrix[0][1] = 1.0;
        matrix[0][2] = 2.0;
        matrix[0][3] = 3.0;
        matrix[1][0] = 4.0;
        matrix[1][1] = 5.0;
        matrix[1][2] = 6.0;
        matrix[1][3] = 7.0;
        matrix[2][0] = 8.0;
        matrix[2][1] = 9.0;
        matrix[2][2] = 10.0;
        matrix[2][3] = 11.0;

        for (int r = 0; r < matrix.length; r++) {
            int indexToRemove = 1;
            System.arraycopy(matrix[r], 0, destMatrix[r], 0, indexToRemove);
            System.arraycopy(matrix[r], indexToRemove + 1, destMatrix[r],
                             indexToRemove,
                             matrix[0].length - (indexToRemove + 1));
        }

        sc.removeVerticalSeam(new int[] { 2, 2, 1, 0, 0, 1, 2 });
        //sc.removeVerticalSeam(new int[] { 1, 1, 1, 1 });

        // Double[][] tMatrixI = sc.castToDouble(sc.transpose(matrix));
        //Integer[][] tMatrixI = Arrays.copyOf(tMatrixO, tMatrixO.length, Integer[].class);

        System.out.println(Arrays.toString(destMatrix));

    }
}
