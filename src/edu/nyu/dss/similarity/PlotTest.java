package edu.nyu.dss.similarity;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.block.BlockBorder;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.Rectangle2D;

public class PlotTest extends JFrame {
    public PlotTest(double[][] data1, double[][] data2) {
        initUI(data1, data2);
    }

    private void initUI(double[][] data1, double[][] data2) {
        XYDataset dataset = createDataset(data1, data2);
        JFreeChart chart = createChart(dataset);
        ChartPanel chartPanel = new ChartPanel(chart);
        add(chartPanel);
        pack();
        setTitle("test title");
        setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
    }

    private XYDataset createDataset(double[][] data1, double[][] data2) {
        XYSeries series1 = new XYSeries("after");
        for(double[] d: data1) {
            series1.add(d[0], d[1]);
        }
        XYSeriesCollection collection = new XYSeriesCollection();
        collection.addSeries(series1);

        XYSeries series2 = new XYSeries("before");
        for(double[] d: data2) {
            series2.add(d[0], d[1]);
        }
        collection.addSeries(series2);

        return collection;
    }

    private JFreeChart createChart(XYDataset dataset) {
        JFreeChart chart = ChartFactory.createXYLineChart(
                "ShapeNet", "indexNode", "radius",
                dataset, PlotOrientation.VERTICAL, true, true, false);
        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
        renderer.setSeriesPaint(0, Color.red);
        renderer.setSeriesStroke(0, new BasicStroke(1.0f));
//        renderer.setSeriesOutlineStroke(0, new BasicStroke(0.1f));
        Shape shape0 = new Rectangle2D.Double(-1, -1, 0.5, 0.5);
        renderer.setSeriesShape(0, shape0);
        renderer.setSeriesPaint(1, Color.blue);
        renderer.setSeriesStroke(1, new BasicStroke(1.0f));
//        renderer.setSeriesOutlineStroke(1, new BasicStroke(0.1f));
        Shape shape1 = new Rectangle2D.Double(-1, -1, 1, 1);
        renderer.setSeriesShape(1, shape1);
        Font font1 = new Font("??", Font.BOLD, 18);
        Font font2 = new Font("??", Font.BOLD, 14);
        XYPlot plot = chart.getXYPlot();
        NumberAxis numberAxis = (NumberAxis) plot.getRangeAxis();
        numberAxis.setTickLabelFont(font2);
        numberAxis.setLabelFont(font1);
        ValueAxis valueAxis = plot.getDomainAxis();
        valueAxis.setTickLabelFont(font2);
        valueAxis.setLabelFont(font1);
        plot.setRenderer(renderer);
        plot.setBackgroundPaint(Color.WHITE);
        plot.setOutlinePaint(Color.black);
        plot.setOutlineStroke(new BasicStroke(1.5f));
        chart.getLegend().setItemFont(font2);
        chart.getLegend().setFrame(BlockBorder.NONE);
        return chart;
    }
}
