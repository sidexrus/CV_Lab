import ui.Gui;

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.geom.AffineTransform;
import java.awt.image.AffineTransformOp;
import java.awt.image.BufferedImage;
import java.awt.image.DataBufferByte;
import java.awt.image.DataBufferInt;
import java.io.File;
import java.util.ArrayList;

public class Main {

    public static void main(String[] args) throws Exception {

        BufferedImage source = ImageIO.read(new File("materials//photo.png"));
        Image src = new Image(ImageIO.read(new File("materials//photo.png")));
        Image img = new Image(ImageIO.read(new File("materials//photo.png")));
        Image img3 = new Image(ImageIO.read(new File("materials//photo.png")));


        //lab1
        //img.DerivativeX(true);
        //BufferedImage bufimg = new BufferedImage(img.getWidth(), img.getHeight(), BufferedImage.TYPE_BYTE_GRAY);
        //bufimg.getRaster().setSamples(0,0, img.getWidth(), img.getHeight(),0, img.Sobel());
        //new Gui().show(img.getImage());
        //new Gui().show(src.getImage());

        //lab2
        ImagePyramid pyr = new ImagePyramid(src);
        pyr.CreatePyramid(5,5,0.5, 1.6);
        pyr.Save_Pyramid();

        //lab3
        InterestPoint ip = new InterestPoint(img3);
        //new Gui().show(ip.HarrisMap(2));
        //new Gui().show(ip.Harris(0.01, 2, 40));
        //new Gui().show(ip.Moravek(0.06, 2, 150));

        //lab4

        Image f = new Image(ImageIO.read(new File("materials//photo.png")));
        InterestPoint ip41 = new InterestPoint(f);
        ip41.Harris(0.05, 2, 30);
        Descriptor[] d1 = DescriptorCreator.getDescriptors(f, ip41.Points, 8, 8, 16);


        Image s = new Image(ImageIO.read(new File("materials//photo.png")));
        InterestPoint ip42 = new InterestPoint(s);
        ip42.Harris(0.05, 2, 30);
        Descriptor[] d2 = DescriptorCreator.getDescriptors(f, ip42.Points, 8, 8, 16);

        BufferedImage res = img3.GlueImages(ip41.res, ip42.res);
        //Image result = new Image(res);

        ArrayList<DescriptorCreator.Vector> similar = DescriptorCreator.findSimilar(d1, d2, (0.8));

        new Gui(res.getWidth(), res.getHeight()).show(Image.DrawLines(similar,f.getWidth(), f.getHeight(), res));


        //lab5

        Image f5 = new Image(ImageIO.read(new File("materials//photo.png")));
        InterestPoint ip51 = new InterestPoint(f5);
        ip51.Harris(0.02, 2, 80);
        Descriptor[] d15 = DescriptorCreator.getDescriptorsInvRotation(f5, ip51.Points, 8, 8, 16);


        BufferedImage r = ImageIO.read(new File("materials//photo.png"));
        r = Image.Rotate(r,45);
        Image s5 = new Image(r);
        InterestPoint ip52 = new InterestPoint(s5);
        ip52.Harris(0.01, 2, 80);
        Descriptor[] d25 = DescriptorCreator.getDescriptorsInvRotation(s5, ip52.Points, 8, 8, 16);

        BufferedImage res5 = img3.GlueImages(ip51.res, ip52.res);
        //Image result = new Image(res);

        ArrayList<DescriptorCreator.Vector> similar5 = DescriptorCreator.findSimilar(d15, d25, (0.8));
        //new Gui(res.getWidth(), res.getHeight()).show(res);
        new Gui(res5.getWidth(), res5.getHeight()).show(Image.DrawLines(similar5,f5.getWidth(), f5.getHeight(), res5));

        /*
        //lab6
        BufferedImage n = ImageIO.read(new File("materials//photo.png"));
        BufferedImage scaledImage = new BufferedImage(n.getWidth()/2, n.getHeight()/2, BufferedImage.TYPE_INT_ARGB);
        Graphics2D graphics2D = scaledImage.createGraphics();
        graphics2D.setRenderingHint(RenderingHints.KEY_INTERPOLATION,RenderingHints.VALUE_INTERPOLATION_BILINEAR);
        graphics2D.drawImage(n, 0, 0, n.getWidth()/2, n.getHeight()/2, null);
        graphics2D.dispose();
        BufferedImage r6 = Image.Rotate(scaledImage, 20);
        Image s6 = new Image(r6);

        new Gui(400, 400).show(r6);
*/

        /*

        new Gui(r.getWidth(), r.getHeight()).show(r);
*/
    }
}
