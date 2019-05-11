import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;

public class ImagePyramid {
    public class Item{
        BufferedImage image;
        double sigma;
        double effective_sigma;
        int size;
        int octave;
        int scale;
        double[] dog_image;
        double[] img;
        public Item(double[] _img, BufferedImage _image, double[] _dog_image, double _sigma, double _effsigma, int _size, int _octave, int _scale){
            image = _image;
            img = _img;
            sigma = _sigma;
            effective_sigma = _effsigma;
            size = _size;
            octave = _octave;
            scale = _scale;
            dog_image = _dog_image;
        }

        public Item(double[] matrix ,BufferedImage _img, double _sigma, double _effsigma, int _size, int _octave, int _scale){
            img = matrix;
            image = _img;
            sigma = _sigma;
            effective_sigma = _effsigma;
            size = _size;
            octave = _octave;
            scale = _scale;
        }
    }
    Image SourceImage;
    ArrayList<Item> Pyramid;
    ArrayList<Item> dogs;


    public ImagePyramid(Image img){
        Pyramid = new ArrayList<Item>();
        dogs = new ArrayList<>();
        SourceImage = img;
    }

    public void L(double sigma){
        SourceImage.GaussianBlur(sigma);
    }

    public void CreatePyramid(int Octave_num, int Scale_num, double sigma0, double sigma1){
        double Sigma_Scale = sigma1;
        double Sigma_Eff = sigma1;
        SourceImage.GaussianBlur(getDeltaSigma(sigma0, sigma1));
        BufferedImage new_img = SourceImage.deepCopy();
        Pyramid.add(new Item(SourceImage.matrix.clone(), new_img, Sigma_Scale, Sigma_Eff,new_img.getHeight()*new_img.getWidth(),
                0, 0));
        double k = Math.pow(2, 1.0 / Scale_num);

        for(int i=0; i<Octave_num; i++){
            for(int j=0; j<Scale_num; j++){
                double sigmaScalePrev = Sigma_Scale;
                Sigma_Scale = sigma1 * Math.pow(k, j + 1);
                double deltaSigma = getDeltaSigma(sigmaScalePrev, Sigma_Scale);
                Sigma_Eff = Sigma_Scale * Math.pow(2, i);
                //SourceImage.GaussianBlur(deltaSigma);
                L(deltaSigma);
                new_img = SourceImage.deepCopy();
                Pyramid.add(new Item(SourceImage.matrix.clone(), new_img, Sigma_Scale, Sigma_Eff, new_img.getHeight()*new_img.getWidth(),
                        i, j+1));
            }
            //SourceImage.Downsampling();
            SourceImage.bilinearHalfReduce();
            new_img = SourceImage.deepCopy();
            Sigma_Scale = sigma1;
            Pyramid.add(new Item(SourceImage.matrix.clone(), new_img, Sigma_Scale, Sigma_Eff, new_img.getHeight()*new_img.getWidth(),
                    i+1, 0));
        }

        /* Constructs DOGs */
        for (int i = 1; i < Pyramid.size(); i++)
        {
            if (Pyramid.get(i - 1).size == Pyramid.get(i).size)
            {
                Item item = Pyramid.get(i - 1);
                //Image img1 = new Image(Pyramid.get(i).img);
                //Image img2 = new Image(item.img);
                double[] img1 = Pyramid.get(i).img;
                double[] img2 = Pyramid.get(i-1).img;
                double[] dog_matrix = new double[img1.length];
                for(int j=0; j<img1.length; j++)
                    dog_matrix[j] = img1[j] - img2[j];
                Item dog = new Item(Pyramid.get(i).img, Pyramid.get(i).image, dog_matrix, item.sigma, item.effective_sigma, dog_matrix.length,
                        item.octave, item.scale);
                dogs.add(dog);
            }
        }
    }

    public double getDeltaSigma(double sigmaPrev, double sigmaNext) {
        return Math.sqrt(sigmaNext * sigmaNext - sigmaPrev * sigmaPrev);
    }

    public void Save_Pyramid(){
        for(int i=0; i<Pyramid.size(); i++)
            try {
                // retrieve image
                String sigma = new DecimalFormat("#0.000").format(Pyramid.get(i).sigma);
                String sigma_eff = new DecimalFormat("#0.000").format(Pyramid.get(i).effective_sigma);
                File outputfile = new File("Lab2_out/" + (i+1) + "_sigma=" + sigma + "_effsigma="+ sigma_eff + ".jpg");
                ImageIO.write(Pyramid.get(i).image, "jpg", outputfile);
            }
            catch (IOException e) {
            }
    }

}
