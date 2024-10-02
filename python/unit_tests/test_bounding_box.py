import unittest
from misc.bounding_box import BoundingBox


class TestBoundingBox(unittest.TestCase):
    def test_constructor(self):
        box = BoundingBox()
        self.assertTrue(all(v is None for v in [box.zmin, box.zmax, box.ymin, box.ymax, box.xmin, box.xmax]))
        zmin, zmax = 1, 2
        ymin, ymax = 3, 4
        xmin, xmax = 5, 6
        box = BoundingBox(zmin=zmin, zmax=zmax, ymin=ymin, ymax=ymax, xmin=xmin, xmax=xmax)
        self.assertEqual(box.zmin, zmin)
        self.assertEqual(box.zmax, zmax)
        self.assertEqual(box.ymin, ymin)
        self.assertEqual(box.ymax, ymax)
        self.assertEqual(box.xmin, xmin)
        self.assertEqual(box.xmax, xmax)
        box_copy = BoundingBox.from_other(box)
        self.assertEqual(box_copy.zmin, zmin)
        self.assertEqual(box_copy.zmax, zmax)
        self.assertEqual(box_copy.ymin, ymin)
        self.assertEqual(box_copy.ymax, ymax)
        self.assertEqual(box_copy.xmin, xmin)
        self.assertEqual(box_copy.xmax, xmax)

    def test_empty(self):
        self.assertFalse(BoundingBox(1, 2, 3, 4, 5, 6).empty())
        self.assertTrue(BoundingBox(None, 2, 3, 4, 5, 6).empty())
        self.assertTrue(BoundingBox().empty())
        self.assertTrue(BoundingBox(1, 2, 3, 3, 5, 6).empty())
        self.assertTrue(BoundingBox(1, 2, 3, 4, 5, -1).empty())

    def test_area(self):
        self.assertEqual(BoundingBox(1, 2, 3, 4).area(), 0)
        self.assertEqual(BoundingBox(2, 1, 3, 4, 5, 5).area(), 0)
        self.assertEqual(BoundingBox(1, 2, 3, 4, 5, 6 ).area(), 1)

    def test_enclose(self):
        box0, box1 = BoundingBox(), BoundingBox(1, 2)
        self.assertTrue(box0.enclose(box1).empty())
        box0, box1 = BoundingBox(1, 3), BoundingBox(1, 2, 3, 4, 5, 6)
        self.assertEqual(box0.enclose(box1), box1)
        self.assertEqual(box1.enclose(box0), box1)
        box0, box1 = BoundingBox(1, 2, 3, 4, 5, 6), BoundingBox(1, 2, 3, 4, 5, 6)
        self.assertEqual(box0.enclose(box1), box0)
        box0, box1 = BoundingBox(1, 2, 3, 5, 3, 8), BoundingBox(-1, 3, 2, 4, 5, 6)
        self.assertEqual(box0.enclose(box1), BoundingBox(-1, 3, 2, 5, 3, 8))
        self.assertEqual(box0, BoundingBox(1, 2, 3, 5, 3, 8))

    def test_intersect(self):
        box0, box1 = BoundingBox(), BoundingBox(1, 2, 3, 4, 5, 6)
        self.assertTrue(box0.intersect(box1).empty())
        self.assertTrue(box1.intersect(box0).empty())
        box0, box1 = BoundingBox(2, 5, 1, 6, -1, 3), BoundingBox(3, 7, 1, 6, 0, 1)
        self.assertEqual(box0.intersect(box1), BoundingBox(3, 5, 1, 6, 0, 1))
        self.assertEqual(box0, BoundingBox(2, 5, 1, 6, -1, 3))
        box0, box1 = BoundingBox(2, 5, 1, 6, -1, 3), BoundingBox(3, 7, 6, 7, 0, 1)
        self.assertTrue(box0.intersect(box1).empty())



