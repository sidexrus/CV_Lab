import javax.imageio.ImageIO;
import java.awt.image.*;
import java.io.ByteArrayOutputStream;
import java.io.Console;
import java.io.IOException;

public class Image {

    private BufferedImage bufimg;
    private int height;
    private int width;
    private int[] matrix;

    public Image(BufferedImage img) throws IOException {
        height = img.getHeight();
        width = img.getWidth();
        matrix = new int[height*width];

        int r, g, b, color;
        for (int i=0; i<height; i++)
            for(int j=0; j<width; j++){
                color = img.getRGB(j, i);
                r = (color & 0xff0000) >> 16;
                g = (color & 0xff00) >> 8;
                b = color & 0xff;
                matrix[i*width+j] = (int)(0.299*r + 0.587*g + 0.114*b);
            }

        bufimg = new BufferedImage(width, height, BufferedImage.TYPE_BYTE_GRAY);
        bufimg.getRaster().setSamples(0,0, width, height,0, matrix);
    }

    public int getWidth() {
        return width;
    }

    public int getHeight() {
        return height;
    }

    public BufferedImage getImage(){
        return  bufimg;
    }

    public BufferedImage deepCopy() {
        ColorModel cm = bufimg.getColorModel();
        boolean isAlphaPremultiplied = cm.isAlphaPremultiplied();
        WritableRaster raster = bufimg.copyData(null);
        return new BufferedImage(cm, raster, isAlphaPremultiplied, null);
    }

    public int[] ImageSupplement (int[] img_matrix, int kernel_w, int kernel_h) {
        int res_width = width + kernel_w - 1;
        int res_height = height + kernel_h - 1;
        int size =  res_height * res_width;
        int[] res_matrix = new int[size];
        for(int i=0;i<size; i++)
            res_matrix[i] = 0;
        int counter = 0;
        for(int i = kernel_h/2; i< res_height - (kernel_h/2); i++)
            for(int j = kernel_w/2; j<res_width - (kernel_w/2);j++)
                res_matrix[i*res_width+j] = img_matrix[counter++];

        for(int i = 0; i<(kernel_h/2); i++)
            for(int j = kernel_w/2; j<res_width - (kernel_w/2); j++){
                res_matrix[i*res_width+j] = res_matrix[kernel_h/2*res_width+j];
                res_matrix[(res_height-i-1)*res_width+j] = res_matrix[(res_height-1-kernel_h/2)*res_width+j];
            }

        for(int i = 0; i<(kernel_w/2); i++)
            for(int j = 0; j<res_height; j++){
                res_matrix[j*res_width+i] = res_matrix[j*res_width+kernel_w/2];
                res_matrix[j*res_width+res_width-i-1] = res_matrix[j*res_width+res_width-kernel_w/2-1];
            }

        return res_matrix;
    }

    public int[] Convolution(double[] kernel, int kernel_w, int kernel_h){
        int[] augmented_img = ImageSupplement(matrix, kernel_w, kernel_h);
        int aug_width = width + kernel_w - 1;
        int aug_height = height + kernel_h - 1;
        int[] res_matrix = new int[width*height];
        int counter = 0;
        int div=0;

        for(int i=0;i<kernel_h*kernel_w;i++)
            div+=kernel[i];

        int kh = kernel_h/2;
        int kw = kernel_w/2;
        for(int i = kh; i < aug_height - kh; i++)
            for(int j = kw; j < aug_width - kw; j++){
                int kernel_counter = 0;
                double value = 0;
                for(int ipix = -kh; ipix <= kh; ipix++)
                    for(int jpix = -kw; jpix <= kw; jpix++) {
                        int position = ipix*aug_width + jpix;
                        value += augmented_img[i * aug_width + j + position]*kernel[kernel_counter++];
                    }
                if(div != 0) value /= div;
                if (value < 0) value = 0;
                if (value > 255) value = 255;
                res_matrix[counter++] = (int)value;
            }
        matrix = res_matrix;
        bufimg.getRaster().setSamples(0,0, width, height,0, matrix);
        return  res_matrix;
    }

    public int[] DerivativeX() {
        double[] kernel = { 1, 0, -1,
                            2, 0, -2,
                            1, 0, -1};
        return Convolution(kernel, 3, 3);
    }

    public int[] DerivativeY() {
        double[] kernel = { 1, 2, 1,
                            0, 0, 0,
                            -1, -2, -1};
        return Convolution(kernel, 3, 3);
    }

    public int[] Sobel(){
        int[] augmented_img = ImageSupplement(matrix, 3, 3);
        int kernel_w = 3;
        int kernel_h = 3;
        int aug_width = width + kernel_w - 1;
        int aug_height = height + kernel_h - 1;
        int[] res_matrix = new int[width*height];
        int counter = 0;

        double[] kernelx = { 1, 0, -1,
                            2, 0, -2,
                            1, 0, -1};

        double[] kernely = { 1, 2, 1,
                            0, 0, 0,
                            -1, -2, -1};

        int kh = kernel_h/2;
        int kw = kernel_w/2;
        for(int i = kh; i < aug_height - kh; i++)
            for(int j = kw; j < aug_width - kw; j++){
                int kernel_counter = 0;
                int valuex = 0;
                int valuey = 0;
                for(int ipix = -kh; ipix <= kh; ipix++)
                    for(int jpix = -kw; jpix <= kw; jpix++) {
                        int position = ipix*aug_width + jpix;
                        valuex += augmented_img[i * aug_width + j + position]*kernelx[kernel_counter];
                        valuey += augmented_img[i * aug_width + j + position]*kernely[kernel_counter++];
                    }
                int value = (int)Math.sqrt(valuex*valuex+valuey*valuey);
                if (value < 0) value = 0;
                if (value > 255) value = 255;
                res_matrix[counter++] = value;
            }
        matrix = res_matrix;
        bufimg.getRaster().setSamples(0,0, width, height,0, matrix);
        return  res_matrix;
    }

    public int[] SeparableFilter(double[] row, double[] column) {
        int length = row.length;
        int size = row.length/2;
        int[] augmented_img = ImageSupplement(matrix, length, length);
        int aug_width = width + length - 1;
        int aug_height = height + length - 1;
        int counter = 0;
        double div=0;
        int[] res_matrix = new int[width*height];

        for(double i:row)
            for(double j:column)
                div+= i*j;

        for(int i = size; i < aug_height - size; i++)
            for(int j = size; j < aug_width - size; j++){
                int kernel_counterY = 0;
                int kernel_counterX = 0;
                double valueX = 0;
                double value = 0;
                //row
                for(int ipix = -size; ipix <= size; ipix++){
                    for(int jpix = -size; jpix <= size; jpix++) {
                        int position = ipix*aug_width + jpix;
                        valueX += augmented_img[i * aug_width + j + position]*row[kernel_counterX++];
                    }
                    //column
                    kernel_counterX = 0;
                    value += valueX*column[kernel_counterY++];
                    valueX = 0;
                }
                if(div != 0) value /= div;
                if (value < 0) value = 0;
                if (value > 255) value = 255;
                res_matrix[counter++] = (int)value;
            }
        return  res_matrix;
    }

    public int[] GaussianBlur(double sigma){
        int size = (int)Math.ceil((double)sigma*2*3);
        if(size%2 == 0)
            size++;
        double[] row = new double[size];
        double[] column = new double[size];
        int s = size/2;
        int counter = 0;
        for(int i=-s; i<=s; i++){
            double f = Math.pow(Math.E, (double)-(i*i)/(2*sigma*sigma));
            row[counter] = f/Math.sqrt(2*Math.PI*sigma);
            column[counter++] = f/Math.sqrt(2*Math.PI*sigma);
        }
        matrix = SeparableFilter(row, column);
        bufimg.getRaster().setSamples(0,0, width, height,0, matrix);

        return SeparableFilter(row, column);
    }

    public int[] Downsampling (){
        int res_width = (int)Math.round((double)width/2);
        int res_height = (int)Math.round((double)height/2);
        int[] result = new int[res_width*res_height];
        int counter = 0;
        for (int i=0; i<height; i+=2)
            for (int j=0; j<width; j+=2)
                result[counter++] = matrix[i*width+j];

        matrix = result;
        width = res_width;
        height = res_height;
        bufimg = new BufferedImage(getWidth(), getHeight(), BufferedImage.TYPE_BYTE_GRAY);
        bufimg.getRaster().setSamples(0,0, res_width, res_height,0, matrix);

        return result;
    }
}
