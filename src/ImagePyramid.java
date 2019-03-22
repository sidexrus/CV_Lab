import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;

public class ImagePyramid {
    public class Item{
        BufferedImage img;
        double sigma;
        double effective_sigma;
        public Item(BufferedImage _img, double _sigma, double _effsigma){
            img = _img;
            sigma = _sigma;
            effective_sigma = _effsigma;
        }
    }
    Image SourceImage;
    ArrayList<Item> Pyramid;

    public ImagePyramid(Image img){
        Pyramid = new ArrayList<>();
        SourceImage = img;
    }

    public void L(double sigma){
        SourceImage.GaussianBlur(sigma);
    }

    public void CreatePyramid(int Octave_num, int Scale_num, double sigma0, double sigma1){
        double Sigma_Scale = sigma1;
        double Sigma_Eff = sigma1;
        SourceImage.GaussianBlur(getDeltaSigma(sigma0, sigma1));
        Pyramid.add(new Item(SourceImage.deepCopy(), Sigma_Scale, Sigma_Eff));
        double k = Math.pow(2, 1.0 / Scale_num);

        for(int i=0; i<Octave_num; i++){
            for(int j=0; j<Scale_num; j++){
                double sigmaScalePrev = Sigma_Scale;
                Sigma_Scale = sigma1 * Math.pow(k, j + 1);
                double deltaSigma = getDeltaSigma(sigmaScalePrev, Sigma_Scale);
                Sigma_Eff = Sigma_Scale * Math.pow(2, i);
                SourceImage.GaussianBlur(deltaSigma);
                Pyramid.add(new Item(SourceImage.deepCopy(), Sigma_Scale, Sigma_Eff));
            }
            SourceImage.Downsampling();
            Pyramid.add(new Item(SourceImage.deepCopy(), Sigma_Scale, Sigma_Eff));
            Sigma_Scale = sigma1;
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
                ImageIO.write(Pyramid.get(i).img, "jpg", outputfile);
            }
            catch (IOException e) {
            }
    }

}
