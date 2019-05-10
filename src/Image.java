import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.geom.AffineTransform;
import java.awt.image.*;
import java.io.ByteArrayOutputStream;
import java.io.Console;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

public class Image {

    public BufferedImage bufimg;
    private int height;
    private int width;
    public double[] matrix;
    public int[] gray;
    static int shift = 30;
    public Image(BufferedImage img) {
        height = img.getHeight();
        width = img.getWidth();
        matrix = new double[height*width];
        gray = new int[height*width];

        int r, g, b, color;
        for (int i=0; i<height; i++)
            for(int j=0; j<width; j++){
                color = img.getRGB(j, i);
                r = (color & 0xff0000) >> 16;
                g = (color & 0xff00) >> 8;
                b = color & 0xff;
                matrix[i*width+j] = (0.299*r + 0.587*g + 0.114*b)/255;
                gray[i*width+j] = (int)(0.299*r + 0.587*g + 0.114*b);

            }
        bufimg = new BufferedImage(width, height, BufferedImage.TYPE_BYTE_GRAY);
        bufimg.getRaster().setSamples(0,0, width, height,0, gray);
    }

    public Image(double[] matr, int height, int width) {
        matrix = matr;
        this.height = height;
        this.width = width;
    }

    public void SetPixel(int x, int y, double value){
        matrix[y*width+x] = value;
    }

    public double GetPixel(int x, int y){
        if (x < 0) x = -x;
        if (y < 0) y = -y;
        if (x >= width) x = 2 * width - x - 1;
        if (y >= height) y = 2 * height - y - 1;
        return matrix[y*width+x];
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

    public double[] ImageSupplement (double[] img_matrix, int kernel_w, int kernel_h) {
        int res_width = width + kernel_w - 1;
        int res_height = height + kernel_h - 1;
        int size =  res_height * res_width;
        double[] res_matrix = new double[size];
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

    public double[] Convolution(double[] kernel, int kernel_w, int kernel_h, boolean norm){
        double[] augmented_img = ImageSupplement(matrix, kernel_w, kernel_h);
        int aug_width = width + kernel_w - 1;
        int aug_height = height + kernel_h - 1;
        double[] res_matrix = new double[width*height];
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
                res_matrix[counter++] = (double)value;
            }
        if(norm) {
            res_matrix = NormalizeTo255(res_matrix);
            bufimg.getRaster().setSamples(0, 0, width, height, 0, res_matrix);
            matrix = NormalizeTo1(res_matrix.clone());
        }
        return  res_matrix;
    }

    public double[] NormalizeTo255(double[] matrix){
        double min = 0;
        double max = 0;
        boolean f = false;

        for(int i=0; i<width*height; i++){
            if(min>matrix[i]) min = matrix[i];
            if(max<matrix[i]) max = matrix[i];
        }

        double DivDifference = 255/(max - min);
        if (min < 0) f = true;
        for(int i=0; i<width*height; i++){
            if (f) matrix[i] += Math.abs(min);
            matrix[i] *=DivDifference;
            if (matrix[i] < 0) matrix[i] = 0;
            if (matrix[i] > 255) matrix[i] = 255;
        }
        //bufimg.getRaster().setSamples(0,0, width, height,0, matrix);
        return  matrix;
    }

    public double[] NormalizeTo1(double[] matrix){
        for(int i=0; i<width*height; i++){
            matrix[i] /=255;
        }
        return  matrix;
    }

    public double[] DerivativeX(boolean norm) {
        double[] kernel = { 1, 0, -1,
                            2, 0, -2,
                            1, 0, -1};
        return Convolution(kernel, 3, 3, norm);
    }

    public double[] DerivativeY(boolean norm) {
        double[] kernel = { 1, 2, 1,
                            0, 0, 0,
                            -1, -2, -1};
        return Convolution(kernel, 3, 3, norm);
    }

    public double[] Sobel(){
        double[] augmented_img = ImageSupplement(matrix, 3, 3);
        int kernel_w = 3;
        int kernel_h = 3;
        int aug_width = width + kernel_w - 1;
        int aug_height = height + kernel_h - 1;
        double[] res_matrix = new double[width*height];
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
                double valuex = 0;
                double valuey = 0;
                for(int ipix = -kh; ipix <= kh; ipix++)
                    for(int jpix = -kw; jpix <= kw; jpix++) {
                        int position = ipix*aug_width + jpix;
                        valuex += augmented_img[i * aug_width + j + position]*kernelx[kernel_counter];
                        valuey += augmented_img[i * aug_width + j + position]*kernely[kernel_counter++];
                    }
                double value = Math.sqrt(valuex*valuex+valuey*valuey);
                //if (value < 0) value = 0;
                //if (value > 255) value = 255;
                res_matrix[counter++] = value;
            }
        res_matrix = NormalizeTo255(res_matrix);
        bufimg.getRaster().setSamples(0,0, width, height,0, res_matrix);
        matrix = NormalizeTo1(res_matrix.clone());
        return  res_matrix;
    }

    public double[] SeparableFilter(double[] row, double[] column) {
        int length = row.length;
        int size = row.length/2;
        double[] augmented_img = ImageSupplement(matrix, length, length);
        int aug_width = width + length - 1;
        int aug_height = height + length - 1;
        int counter = 0;
        double div=0;
        double[] res_matrix = new double[width*height];

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
                //if (value < 0) value = 0;
                //if (value > 255) value = 255;
                res_matrix[counter++] = (double) value;
            }
        res_matrix = NormalizeTo255(res_matrix);
        bufimg.getRaster().setSamples(0,0, width, height,0, res_matrix);
        matrix = NormalizeTo1(res_matrix.clone());
        return  res_matrix;
    }

    public double[] GaussianBlur(double sigma){
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
        //matrix = SeparableFilter(row, column);
        //bufimg.getRaster().setSamples(0,0, width, height,0, matrix);

        return SeparableFilter(row, column);
    }

    public double[] GaussKernel(double sigma){
        int size = (int)Math.ceil((double)sigma*2*3);
        if(size%2 == 0)
            size++;
        double[] matrix = new double[size*size];
        int s = size/2;
        int k=0;
        double sum=0;
        for(int x=-s; x<=s; x++){
            for(int y=-s; y<=s; y++) {
                double f = Math.pow(Math.E, (double)-(x*x+y*y)/(2*sigma*sigma));
                matrix[k] = f/(2*Math.PI*sigma*sigma);
                sum+=matrix[k];
                k++;
            }
        }
        for(int i=0; i<size*size; i++)
            matrix[i]/=sum;

        return  matrix;
    }

    public double[] Downsampling (){
        int res_width = (int)Math.round((double)width/2);
        int res_height = (int)Math.round((double)height/2);
        double[] result = new double[res_width*res_height];
        int counter = 0;
        for (int i=0; i<height; i+=2){
            for (int j=0; j<width; j+=2)
                result[counter++] = matrix[i*width+j];
        }

        width = res_width;
        height = res_height;
        result = NormalizeTo255(result);
        bufimg = new BufferedImage(getWidth(), getHeight(), BufferedImage.TYPE_BYTE_GRAY);
        bufimg.getRaster().setSamples(0,0, res_width, res_height,0, result);

        matrix = NormalizeTo1(result.clone());

        return result;
    }

    public BufferedImage GlueImages(Image img1, Image img2){
        int Width = img1.getWidth() + img2.getWidth() + shift;
        int Height;
        if ((img1.getHeight() > (img2.getHeight() + shift))) {
            Height = img1.getHeight();
        } else {
            Height = img2.getHeight() + shift;
        }
        BufferedImage res = new BufferedImage(Width, Height, BufferedImage.TYPE_BYTE_GRAY);
        width = img1.getWidth();
        height = img1.getHeight();
        res.getRaster().setSamples(0,0, img1.getWidth(), img1.getHeight(),0, NormalizeTo255(img1.matrix));
        width = img2.getWidth();
        height = img2.getHeight();
        res.getRaster().setSamples(img1.getWidth()+shift,shift, img2.getWidth(), img2.getHeight(),0, NormalizeTo255(img2.matrix));

        return res;
    }

    public static BufferedImage DrawLines(ArrayList<DescriptorCreator.Vector> similar, int width, int height, BufferedImage img){
        BufferedImage res = new BufferedImage(img.getWidth(), img.getHeight(), BufferedImage.TYPE_INT_RGB);

        //res.getRaster().setSamples(0,0, img.getWidth(), img.getHeight(),0, img.getRGB(0, 0, img.getWidth(), img.getHeight(), null, 0, img.getWidth()));
        Graphics2D g = res.createGraphics();
        g.drawImage(img, null, 0, 0);
        BasicStroke bs =new BasicStroke(1);
        g.setStroke(bs);
        for(int i=0; i< similar.size(); i++){
            int x1 = similar.get(i).first.getInterPoint().x;
            int y1 = similar.get(i).first.getInterPoint().y;
            int x2 = similar.get(i).second.getInterPoint().x + shift + width;
            int y2 = similar.get(i).second.getInterPoint().y + shift;
            Random rand = new Random();
            float r = rand.nextFloat();
            float gr = rand.nextFloat();
            float b = rand.nextFloat();
            java.awt.Color cc = new java.awt.Color(r, gr, b);
            g.setPaint(cc);
            g.drawLine(x1, y1, x2, y2);
        }
        return res;
    }

    public static BufferedImage Rotate(BufferedImage img, int angle){

        AffineTransform tx = new AffineTransform();
        double theta = angle*Math.PI/180;
        tx.translate((float)img.getWidth() / 2, (float)img.getHeight() / 2);
        tx.rotate(Math.toRadians(angle));
        tx.translate(-(float)img.getWidth() / 2, -(float)img.getHeight() / 2);
        AffineTransformOp op = new AffineTransformOp(tx,
                AffineTransformOp.TYPE_BILINEAR);
        img = op.filter(img, null);

        return img;
    }
}
